*
*     TEST SUBROUTINES FOR NONSMOOTH OPTIMIZATION
*
* SUBROUTINE TIUD06             ALL SYSTEMS                 99/12/01
C PORTABILITY : ALL SYSTEMS
C 90/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
*  INITIATION OF VARIABLES FOR NONLINEAR MINIMAX APPROXIMATION.
*  UNCONSTRAINED DENSE VERSION.
*
* PARAMETERS :
*  IO  N  NUMBER OF VARIABLES.
*  IO  NA  NUMBER OF PARTIAL FUNCTIONS.
*  RO  X(N)  VECTOR OF VARIABLES.
*  RO  FMIN  LOWER BOUND FOR VALUE OF THE OBJECTIVE FUNCTION.
*  RO  XMAX  MAXIMUM STEPSIZE.
*  IO  NEXT  NUMBER OF THE TEST PROBLEM.
*  IO  IEXT  TYPE OF OBJECTIVE FUNCTION. IEXT<0-MAXIMUM OF VALUES.
*         IEXT=0-MAXIMUM OF ABSOLUTE VALUES.
*  IO  IERR  ERROR INDICATOR.
*
      SUBROUTINE TIUD06(N,NA,X,FMIN,XMAX,NEXT,IEXT,IERR)
C     .. Scalar Arguments ..
      DOUBLE PRECISION FMIN,XMAX
      INTEGER IERR,IEXT,N,NA,NEXT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION X(N)
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION Y(123)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION T
      INTEGER I
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,COS,DBLE,EXP,SIN,SQRT
C     ..
C     .. Common blocks ..
      COMMON /EMPR06/Y
C     ..
      FMIN = -1.0D60
      XMAX = 1.0D3
      IEXT = -1
      IERR = 0
      GO TO (10,20,30,40,50,70,90,120,130,
     +       140,150,180,200,220,240,260,280,290,
     +       350,360,370,380,400,420,440) NEXT

   10 IF (N.GE.2 .AND. NA.GE.3) THEN
          N = 2
          NA = 3
          X(1) = 2.0D0
          X(2) = 2.0D0

      ELSE
          IERR = 1
      END IF

      RETURN

   20 IF (N.GE.2 .AND. NA.GE.3) THEN
          N = 2
          NA = 3
          X(1) = 3.0D0
          X(2) = 1.0D0

      ELSE
          IERR = 1
      END IF

      RETURN

   30 IF (N.GE.2 .AND. NA.GE.2) THEN
          N = 2
          NA = 2
          X(1) = 1.41831D0
          X(2) = -4.79462D0
          XMAX = 1.0D0

      ELSE
          IERR = 1
      END IF

      RETURN

   40 IF (N.GE.3 .AND. NA.GE.6) THEN
          N = 3
          NA = 6
          X(1) = 1.0D0
          X(2) = 1.0D0
          X(3) = 1.0D0

      ELSE
          IERR = 1
      END IF

      RETURN

   50 IF (N.GE.4 .AND. NA.GE.4) THEN
          N = 4
          NA = 4
          DO 60 I = 1,N
              X(I) = 0.0D0
   60     CONTINUE

      ELSE
          IERR = 1
      END IF

      RETURN

   70 IF (N.GE.4 .AND. NA.GE.4) THEN
          N = 4
          NA = 4
          DO 80 I = 1,N
              X(I) = 0.0D0
   80     CONTINUE

      ELSE
          IERR = 1
      END IF

      RETURN

   90 IF (N.GE.3 .AND. NA.GE.21) THEN
          N = 3
          NA = 21
          DO 100 I = 1,N
              X(I) = 1.0D0
  100     CONTINUE
          DO 110 I = 1,NA
              T = 1.0D1*DBLE(I-1)/DBLE(NA-1)
              Y(I) = T
              Y(NA+I) = (3.0D0/2.0D1)*EXP(-T) +
     +                  (1.0D0/5.2D1)*EXP(-5.0D0*T) -
     +                  (1.0D0/6.5D1)*EXP(-2.0D0*T)*
     +                  (3.0D0*SIN(2.0D0*T)+1.1D1*COS(2.0D0*T))
  110     CONTINUE
          IEXT = 0

      ELSE
          IERR = 1
      END IF

      RETURN

  120 IF (N.GE.3 .AND. NA.GE.15) THEN
          N = 3
          NA = 15
          X(1) = 1.0D0
          X(2) = 1.0D0
          X(3) = 1.0D0
          Y(1) = 0.14D0
          Y(2) = 0.18D0
          Y(3) = 0.22D0
          Y(4) = 0.25D0
          Y(5) = 0.29D0
          Y(6) = 0.32D0
          Y(7) = 0.35D0
          Y(8) = 0.39D0
          Y(9) = 0.37D0
          Y(10) = 0.58D0
          Y(11) = 0.73D0
          Y(12) = 0.96D0
          Y(13) = 1.34D0
          Y(14) = 2.10D0
          Y(15) = 4.39D0
          IEXT = 0

      ELSE
          IERR = 1
      END IF

      RETURN

  130 IF (N.GE.4 .AND. NA.GE.11) THEN
          N = 4
          NA = 11
          X(1) = 0.25D0
          X(2) = 0.39D0
          X(3) = 4.15D-1
          X(4) = 0.39D0
          Y(1) = 0.1957D0
          Y(2) = 0.1947D0
          Y(3) = 0.1735D0
          Y(4) = 0.1600D0
          Y(5) = 0.0844D0
          Y(6) = 0.0627D0
          Y(7) = 0.0456D0
          Y(8) = 0.0342D0
          Y(9) = 0.0323D0
          Y(10) = 0.0235D0
          Y(11) = 0.0246D0
          Y(12) = 4.0000D0
          Y(13) = 2.0000D0
          Y(14) = 1.0000D0
          Y(15) = 0.5000D0
          Y(16) = 0.2500D0
          Y(17) = 0.1670D0
          Y(18) = 0.1250D0
          Y(19) = 0.1000D0
          Y(20) = 0.0833D0
          Y(21) = 0.0714D0
          Y(22) = 0.0625D0
          IEXT = 0

      ELSE
          IERR = 1
      END IF

      RETURN

  140 IF (N.GE.4 .AND. NA.GE.20) THEN
          N = 4
          NA = 20
          X(1) = 2.5D1
          X(2) = 5.0D0
          X(3) = -5.0D0
          X(4) = -1.0D0
          IEXT = 0

      ELSE
          IERR = 1
      END IF

      RETURN

  150 IF (N.GE.4 .AND. NA.GE.21) THEN
          N = 4
          NA = 21
          DO 160 I = 1,N
              X(I) = 1.0D0
  160     CONTINUE
          DO 170 I = 1,NA
              Y(I) = 0.25D0 + 0.75D0*DBLE(I-1)/DBLE(NA-1)
              Y(NA+I) = SQRT(Y(I))
  170     CONTINUE
          IEXT = 0

      ELSE
          IERR = 1
      END IF

      RETURN

  180 IF (N.GE.4 .AND. NA.GE.21) THEN
          N = 4
          NA = 21
          X(1) = 1.0D0
          X(2) = 1.0D0
          X(3) = -3.0D0
          X(4) = -1.0D0
          DO 190 I = 1,NA
              Y(I) = -0.5D0 + DBLE(I-1)/DBLE(NA-1)
              Y(NA+I) = 1.0D0/ (1.0D0+Y(I))
  190     CONTINUE
          XMAX = 1.0D-1
          IEXT = 0

      ELSE
          IERR = 1
      END IF

      RETURN

  200 IF (N.GE.4 .AND. NA.GE.61) THEN
          N = 4
          NA = 61
          IEXT = 0
          DO 210 I = 1,N
              X(I) = 1.0D0
  210     CONTINUE
          X(3) = 1.0D1
          Y(1) = 1.0D0
          Y(2) = 1.01D0
          Y(3) = 1.02D0
          Y(4) = 1.03D0
          Y(5) = 1.05D0
          Y(6) = 1.075D0
          Y(7) = 1.1D0
          Y(8) = 1.125D0
          Y(9) = 1.15D0
          Y(10) = 1.2D0
          Y(11) = 1.25D0
          Y(12) = 1.3D0
          Y(13) = 1.35D0
          Y(14) = 1.4D0
          Y(15) = 1.5D0
          Y(16) = 1.6D0
          Y(17) = 1.7D0
          Y(18) = 1.8D0
          Y(19) = 1.9D0
          Y(20) = 2.0D0
          Y(21) = 2.1D0
          Y(22) = 2.2D0
          Y(23) = 2.3D0
          Y(24) = 2.5D0
          Y(25) = 2.75D0
          Y(26) = 3.0D0
          Y(27) = 3.25D0
          Y(28) = 3.5D0
          Y(29) = 4.0D0
          Y(30) = 4.5D0
          Y(31) = 5.0D0
          Y(32) = 5.5D0
          Y(33) = 6.0D0
          Y(34) = 6.5D0
          Y(35) = 7.0D0
          Y(36) = 7.5D0
          Y(37) = 8.0D0
          Y(38) = 8.5D0
          Y(39) = 9.0D0
          Y(40) = 10.0D0
          Y(41) = 11.0D0
          Y(42) = 12.0D0
          Y(43) = 13.0D0
          Y(44) = 15.0D0
          Y(45) = 17.5D0
          Y(46) = 20.0D0
          Y(47) = 22.5D0
          Y(48) = 25.0D0
          Y(49) = 30.0D0
          Y(50) = 35.0D0
          Y(51) = 40.0D0
          Y(52) = 50.0D0
          Y(53) = 60.0D0
          Y(54) = 70.0D0
          Y(55) = 80.0D0
          Y(56) = 100.0D0
          Y(57) = 150.0D0
          Y(58) = 200.0D0
          Y(59) = 300.0D0
          Y(60) = 500.0D0
          Y(61) = 1.0D5
          Y(61+1) = 0.97386702052733792831D0
          Y(61+2) = 0.97390711665677071911D0
          Y(61+3) = 0.97394794566286525039D0
          Y(61+4) = 0.97398947529386626621D0
          Y(61+5) = 0.97407451325974368215D0
          Y(61+6) = 0.97418422166965892644D0
          Y(61+7) = 0.97429732692565188272D0
          Y(61+8) = 0.97441344289222034304D0
          Y(61+9) = 0.97453221704823108216D0
          Y(61+10) = 0.97477647977277153145D0
          Y(61+11) = 0.97502785781178233026D0
          Y(61+12) = 0.97528446418205610067D0
          Y(61+13) = 0.97554472005909873148D0
          Y(61+14) = 0.97580730389916439626D0
          Y(61+15) = 0.97633521198091785788D0
          Y(61+16) = 0.97686134356195586299D0
          Y(61+17) = 0.97738094095418268249D0
          Y(61+18) = 0.97789073928751194169D0
          Y(61+19) = 0.97838854811088140808D0
          Y(61+20) = 0.97887295363155439576D0
          Y(61+21) = 0.97934310478576951385D0
          Y(61+22) = 0.97979855827226762515D0
          Y(61+23) = 0.98023916551033862691D0
          Y(61+24) = 0.98107624468416045728D0
          Y(61+25) = 0.98204290774765289406D0
          Y(61+26) = 0.98292719363632655668D0
          Y(61+27) = 0.98373656564197279264D0
          Y(61+28) = 0.98447846610682328991D0
          Y(61+29) = 0.98578713114264981186D0
          Y(61+30) = 0.98690124654380846379D0
          Y(61+31) = 0.98785879054855173380D0
          Y(61+32) = 0.98868928566806726978D0
          Y(61+33) = 0.98941568049711884384D0
          Y(61+34) = 0.99005592865089067038D0
          Y(61+35) = 0.99062420259214811899D0
          Y(61+36) = 0.99113180018738487730D0
          Y(61+37) = 0.99158781685339306121D0
          Y(61+38) = 0.99199964493176098231D0
          Y(61+39) = 0.99237334707422899195D0
          Y(61+40) = 0.99302559755582945576D0
          Y(61+41) = 0.99357562712206729735D0
          Y(61+42) = 0.99404560031581354300D0
          Y(61+43) = 0.99445173790980305195D0
          Y(61+44) = 0.99511816085114882367D0
          Y(61+45) = 0.99575584307408838284D0
          Y(61+46) = 0.99624640327264396775D0
          Y(61+47) = 0.99663543022201287399D0
          Y(61+48) = 0.99695146031888813172D0
          Y(61+49) = 0.99743367936799001685D0
          Y(61+50) = 0.99778424120023198554D0
          Y(61+51) = 0.99805056960591223604D0
          Y(61+52) = 0.99842841443786596919D0
          Y(61+53) = 0.99868358857261655169D0
          Y(61+54) = 0.99886748198687248566D0
          Y(61+55) = 0.99900629944600342584D0
          Y(61+56) = 0.99920194660435455419D0
          Y(61+57) = 0.99946519560889341627D0
          Y(61+58) = 0.99959785208794891934D0
          Y(61+59) = 0.99973120214935885075D0
          Y(61+60) = 0.99983838442420395745D0
          Y(61+61) = 0.999999189398046846077D0

      ELSE
          IERR = 1
      END IF

      RETURN

  220 IF (N.GE.5 .AND. NA.GE.21) THEN
          N = 5
          NA = 21
          DO 230 I = 1,N
              X(I) = 0.0D0
  230     CONTINUE
          X(1) = 0.5D0
          IEXT = 0

      ELSE
          IERR = 1
      END IF

      RETURN

  240 IF (N.GE.5 .AND. NA.GE.30) THEN
          N = 5
          NA = 30
          X(1) = 0.0D0
          X(2) = -1.0D0
          X(3) = 1.0D1
          X(4) = 1.0D0
          X(5) = 1.0D1
          DO 250 I = 1,NA
              Y(I) = -1.0D0 + 2.0D0*DBLE(I-1)/DBLE(NA-1)
              T = 8.0D0*Y(I)
              Y(NA+I) = SQRT((T-1.0D0)**2+1.0D0)*ATAN(T)/T
  250     CONTINUE
          IEXT = 0

      ELSE
          IERR = 1
      END IF

      RETURN

  260 IF (N.GE.6 .AND. NA.GE.51) THEN
          N = 6
          NA = 51
          X(1) = 2.0D0
          X(2) = 2.0D0
          X(3) = 7.0D0
          X(4) = 0.0D0
          X(5) = -2.0D0
          X(6) = 1.0D0
          DO 270 I = 1,NA
              T = 0.1D0*DBLE(I-1)
              Y(I) = 0.5D0*EXP(-T) - EXP(-2.0D0*T) +
     +               0.5D0*EXP(-3.0D0*T) + 1.5D0*EXP(-1.5D0*T)*
     +               SIN(7.0D0*T) + EXP(-2.5D0*T)*SIN(5.0D0*T)
  270     CONTINUE
          IEXT = 0

      ELSE
          IERR = 1
      END IF

      RETURN

  280 IF (N.GE.6 .AND. NA.GE.11) THEN
          N = 6
          NA = 11
          X(1) = 0.8D0
          X(2) = 1.5D0
          X(3) = 1.2D0
          X(4) = 3.0D0
          X(5) = 0.8D0
          X(6) = 6.0D0
          Y(1) = 0.5D0
          Y(2) = 0.6D0
          Y(3) = 0.7D0
          Y(4) = 0.77D0
          Y(5) = 0.9D0
          Y(6) = 1.0D0
          Y(7) = 1.1D0
          Y(8) = 1.23D0
          Y(9) = 1.3D0
          Y(10) = 1.4D0
          Y(11) = 1.5D0

      ELSE
          IERR = 1
      END IF

      RETURN

  290 IF (N.GE.9 .AND. NA.GE.41) THEN
          N = 9
          NA = 41
          X(1) = 0D0
          X(2) = 1D0
          X(3) = 0D0
          X(4) = -1.5D-1
          X(5) = 0D0
          X(6) = -6.8D-1
          X(7) = 0D0
          X(8) = -7.2D-1
          X(9) = 3.7D-1
          DO 300 I = 1,6
              Y(I) = 1D-2* (I-1)
  300     CONTINUE
          DO 310 I = 7,20
              Y(I) = 3D-2* (I-7) + 7D-2
  310     CONTINUE
          Y(21) = 5D-1
          DO 320 I = 22,35
              Y(I) = 3D-2* (I-22) + 54D-2
  320     CONTINUE
          DO 330 I = 36,41
              Y(I) = 1D-2* (I-36) + 95D-2
  330     CONTINUE
          DO 340 I = 1,41
              Y(41+I) = COS(Y(I)*3.14159265358979324D0)
              Y(82+I) = SIN(Y(I)*3.14159265358979324D0)
  340     CONTINUE
          IEXT = 0

      ELSE
          IERR = 1
      END IF

      RETURN

  350 IF (N.GE.7 .AND. NA.GE.5) THEN
          N = 7
          NA = 5
          X(1) = 1.0D0
          X(2) = 2.0D0
          X(3) = 0.0D0
          X(4) = 4.0D0
          X(5) = 0.0D0
          X(6) = 1.0D0
          X(7) = 1.0D0

      ELSE
          IERR = 1
      END IF

      RETURN

  360 IF (N.GE.10 .AND. NA.GE.9) THEN
          N = 10
          NA = 9
          X(1) = 2.0D0
          X(2) = 3.0D0
          X(3) = 5.0D0
          X(4) = 5.0D0
          X(5) = 1.0D0
          X(6) = 2.0D0
          X(7) = 7.0D0
          X(8) = 3.0D0
          X(9) = 6.0D0
          X(10) = 1.0D1

      ELSE
          IERR = 1
      END IF

      RETURN

  370 IF (N.GE.20 .AND. NA.GE.18) THEN
          N = 20
          NA = 18
          X(1) = 2.0D0
          X(2) = 3.0D0
          X(3) = 5.0D0
          X(4) = 5.0D0
          X(5) = 1.0D0
          X(6) = 2.0D0
          X(7) = 7.0D0
          X(8) = 3.0D0
          X(9) = 6.0D0
          X(10) = 1.0D1
          X(11) = 2.0D0
          X(12) = 2.0D0
          X(13) = 6.0D0
          X(14) = 1.5D1
          X(15) = 1.0D0
          X(16) = 2.0D0
          X(17) = 1.0D0
          X(18) = 2.0D0
          X(19) = 1.0D0
          X(20) = 3.0D0

      ELSE
          IERR = 1
      END IF

      RETURN

  380 IF (N.GE.10 .AND. NA.GE.2) THEN
          N = 10
          NA = 2
          DO 390 I = 1,N
              X(I) = 0.1D0
  390     CONTINUE
          X(1) = 1.0D2
          XMAX = 1.0D1

      ELSE
          IERR = 1
      END IF

      RETURN

  400 IF (N.GE.11 .AND. NA.GE.10) THEN
          N = 11
          NA = 10
          DO 410 I = 1,N
              X(I) = 1.0D0
  410     CONTINUE

      ELSE
          IERR = 1
      END IF

      RETURN

  420 IF (N.GE.20 .AND. NA.GE.31) THEN
          N = 20
          NA = 31
          DO 430 I = 1,N
              X(I) = 0.0D0
  430     CONTINUE
          IEXT = 0

      ELSE
          IERR = 1
      END IF

      RETURN

  440 IF (N.GE.11 .AND. NA.GE.65) THEN
          N = 11
          NA = 65
          X(1) = 1.3D0
          X(2) = 6.5D-1
          X(3) = 6.5D-1
          X(4) = 0.7D0
          X(5) = 0.6D0
          X(6) = 3.0D0
          X(7) = 5.0D0
          X(8) = 7.0D0
          X(9) = 2.0D0
          X(10) = 4.5D0
          X(11) = 5.5D0
          Y(1) = 1.366D0
          Y(2) = 1.191D0
          Y(3) = 1.112D0
          Y(4) = 1.013D0
          Y(5) = 0.991D0
          Y(6) = 0.885D0
          Y(7) = 0.831D0
          Y(8) = 0.847D0
          Y(9) = 0.786D0
          Y(10) = 0.725D0
          Y(11) = 0.746D0
          Y(12) = 0.679D0
          Y(13) = 0.608D0
          Y(14) = 0.655D0
          Y(15) = 0.616D0
          Y(16) = 0.606D0
          Y(17) = 0.602D0
          Y(18) = 0.626D0
          Y(19) = 0.651D0
          Y(20) = 0.724D0
          Y(21) = 0.649D0
          Y(22) = 0.649D0
          Y(23) = 0.694D0
          Y(24) = 0.644D0
          Y(25) = 0.624D0
          Y(26) = 0.661D0
          Y(27) = 0.612D0
          Y(28) = 0.558D0
          Y(29) = 0.553D0
          Y(30) = 0.495D0
          Y(31) = 0.500D0
          Y(32) = 0.423D0
          Y(33) = 0.395D0
          Y(34) = 0.375D0
          Y(35) = 0.372D0
          Y(36) = 0.391D0
          Y(37) = 0.396D0
          Y(38) = 0.405D0
          Y(39) = 0.428D0
          Y(40) = 0.429D0
          Y(41) = 0.523D0
          Y(42) = 0.562D0
          Y(43) = 0.607D0
          Y(44) = 0.653D0
          Y(45) = 0.672D0
          Y(46) = 0.708D0
          Y(47) = 0.633D0
          Y(48) = 0.668D0
          Y(49) = 0.645D0
          Y(50) = 0.632D0
          Y(51) = 0.591D0
          Y(52) = 0.559D0
          Y(53) = 0.597D0
          Y(54) = 0.625D0
          Y(55) = 0.739D0
          Y(56) = 0.710D0
          Y(57) = 0.729D0
          Y(58) = 0.720D0
          Y(59) = 0.636D0
          Y(60) = 0.581D0
          Y(61) = 0.428D0
          Y(62) = 0.292D0
          Y(63) = 0.162D0
          Y(64) = 0.098D0
          Y(65) = 0.054D0
          XMAX = 1.0D1
          IEXT = 0

      ELSE
          IERR = 1
      END IF

      RETURN

      END
* SUBROUTINE TAFU06             ALL SYSTEMS                 99/12/01
C PORTABILITY : ALL SYSTEMS
C 91/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
*  VALUES OF PARTIAL FUNCTIONS IN THE MINIMAX CRITERION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  II  KA  INDEX OF THE PARTIAL FUNCTION.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  FA  VALUE OF THE PARTIAL FUNCTION AT THE
*          SELECTED POINT.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*
      SUBROUTINE TAFU06(N,KA,X,FA,NEXT)
C     .. Parameters ..
      DOUBLE PRECISION PI
      PARAMETER (PI=3.14159265358979323846D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION FA
      INTEGER KA,N,NEXT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION X(N)
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION Y(123)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX C1,C2,C3
      DOUBLE PRECISION BETA,T,X1,X2,X3,X4,X5,X6,X7,X8
      INTEGER I,J
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX CA(4),CB(4)
      DOUBLE PRECISION XA(3),XB(3),XC(3)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,CDABS,CMPLX,COS,DBLE,EXP,MIN,SIN,SQRT
C     ..
C     .. Common blocks ..
      COMMON /EMPR06/Y
C     ..
      GO TO (10,50,60,70,140,190,240,250,260,
     +       270,280,290,300,310,330,340,350,380,
     +       390,450,550,660,680,700,720) NEXT

   10 X1 = X(1)*X(1)
      X2 = X(2)*X(2)
      X3 = X(1) + X(1)
      X4 = X(2) + X(2)
      GO TO (20,30,40) KA

   20 FA = X1 + X2*X2
      RETURN

   30 FA = 8.0D0 - 4.0D0* (X(1)+X(2)) + X1 + X2
      RETURN

   40 FA = 2.0D0*EXP(X(2)-X(1))
      RETURN

   50 X1 = 1.0D1*X(1)/ (X(1)+1.0D-1)
      X2 = 2.0D0*X(2)**2
      IF (KA.EQ.1) THEN
          FA = 0.5D0* (X(1)+X1+X2)

      ELSE IF (KA.EQ.2) THEN
          FA = 0.5D0* (-X(1)+X1+X2)

      ELSE IF (KA.EQ.3) THEN
          FA = 0.5D0* (X(1)-X1+X2)
      END IF

      RETURN

   60 X1 = X(1)**2 + X(2)**2
      X2 = SQRT(X1)
      IF (KA.EQ.1) THEN
          FA = (X(1)-X2*COS(X2))**2 + 5.0D-3*X1

      ELSE IF (KA.EQ.2) THEN
          FA = (X(2)-X2*SIN(X2))**2 + 5.0D-3*X1
      END IF

      RETURN

   70 GO TO (80,90,100,110,120,130) KA

   80 FA = X(1)**2 + X(2)**2 + X(3)**2 - 1.0D0
      RETURN

   90 FA = X(1)**2 + X(2)**2 + (X(3)-2.0D0)**2
      RETURN

  100 FA = X(1) + X(2) + X(3) - 1.0D0
      RETURN

  110 FA = X(1) + X(2) - X(3) + 1.0D0
      RETURN

  120 FA = 2.0D0* (X(1)**3+3.0D0*X(2)**2+ (5.0D0*X(3)-X(1)+1.0D0)**2)
      RETURN

  130 FA = X(1)**2 - 9.0D0*X(3)
      RETURN

  140 X1 = X(1)*X(1)
      X2 = X(2)*X(2)
      X3 = X(3)*X(3)
      X4 = X(4)*X(4)
      X5 = X(1) + X(1)
      X6 = X(2) + X(2)
      X7 = X(3) + X(3)
      X8 = X(4) + X(4)
      FA = X1 + X2 + X3 + X3 + X4 - 5.0D0* (X(1)+X(2)) - 2.1D1*X(3) +
     +     7.0D0*X(4)
  150 GO TO (320,160,170,180) KA

  160 FA = FA + 1.0D1* (X1+X2+X3+X4+X(1)-X(2)+X(3)-X(4)-8.0D0)
      RETURN

  170 FA = FA + 1.0D1* (X1+X2+X2+X3+X4+X4-X(1)-X(4)-1.0D1)
      RETURN

  180 FA = FA + 1.0D1* (X1+X2+X3+X5-X(2)-X(4)-5.0D0)
      RETURN

  190 X1 = X(1) - (X(4)+1.0D0)**4
      X2 = X1*X1
      X3 = X(2) - X2*X2
      X4 = X3*X3
      FA = X2 + X4 + 2.0D0*X(3)**2 + X(4)**2 - 5.0D0* (X1+X3) -
     +     2.1D1*X(3) + 7.0D0*X(4)
      GO TO (200,210,220,230) KA

  200 CONTINUE
      RETURN

  210 FA = FA + 1.0D1* (X2+X4+X(3)**2+X(4)**2+X1-X3+X(3)-X(4)-8.0D0)
      RETURN

  220 FA = FA + 1.0D1* (X2+2.0D0*X4+X(3)**2+2.0D0*X(4)**2-X1-X(4)-1.0D1)
      RETURN

  230 FA = FA + 1.0D1* (X2+X4+X(3)**2+2.0D0*X1-X3-X(4)-5.0D0)
      RETURN

  240 T = Y(KA)
      FA = (X(3)/X(2))*EXP(-X(1)*T)*SIN(X(2)*T) - Y(KA+21)
      RETURN

  250 FA = Y(KA) - X(1) - DBLE(KA)/ (DBLE(16-KA)*X(2)+
     +     DBLE(MIN(KA,16-KA))*X(3))
      RETURN

  260 T = Y(KA+11)
      FA = Y(KA) - X(1)*T* (T+X(2))/ ((T+X(3))*T+X(4))
      RETURN

  270 T = 0.2D0*DBLE(KA)
      FA = (X(1)+X(2)*T-EXP(T))**2 + (X(3)+X(4)*SIN(T)-COS(T))**2
      RETURN

  280 T = Y(KA)
      FA = X(4) - ((X(1)*T+X(2))*T+X(3))**2 - Y(KA+21)
      RETURN

  290 T = Y(KA)
      FA = X(1)*EXP(X(3)*T) + X(2)*EXP(X(4)*T) - Y(KA+21)
      RETURN

  300 T = Y(KA)
      FA = X(1)*ABS((T+X(2)+1.0D0/ (X(3)*T+X(4)))/
     +     ((T+1.0D0)*Y(61+KA)))** (T+5.0D-1) - 1.0D0
      RETURN

  310 T = 0.1D0*DBLE(KA-1) - 1.0D0
      X1 = X(1) + T*X(2)
      X2 = 1.0D0/ (1.0D0+T* (X(3)+T* (X(4)+T*X(5))))
      X3 = X1*X2 - EXP(T)
      FA = X3
  320 RETURN

  330 T = Y(KA)
      FA = (X(1)+T* (X(2)+T*X(3)))/ (1.0D0+T* (X(4)+T*X(5))) - Y(KA+30)
      RETURN

  340 T = 0.1D0*DBLE(KA-1)
      FA = X(1)*EXP(-X(2)*T)*COS(X(3)*T+X(4)) + X(5)*EXP(-X(6)*T) -
     +     Y(KA)
      RETURN

  350 BETA = 0.5D0*PI*Y(KA)
      DO 360 I = 1,3
          J = I + I
          XA(I) = X(J-1)
          XB(I) = X(J)
  360 CONTINUE
      CA(4) = CMPLX(1.0D0,0.0D0)
      CB(4) = 1.0D1*CA(4)
      DO 370 J = 1,3
          I = 4 - J
          XC(I) = BETA*XA(I)
          T = XC(I)
          X1 = COS(T)
          X2 = SIN(T)
          C1 = CMPLX(X1,0.0D0)
          C2 = CMPLX(0.0D0, (X2*XB(I)))
          C3 = CMPLX(0.0D0, (X2/XB(I)))
          CB(I) = C1*CB(I+1) + C2*CA(I+1)
          CA(I) = C3*CB(I+1) + C1*CA(I+1)
  370 CONTINUE
      C1 = -CA(1)
      C2 = CB(1) - C1
      C3 = 1.0D0 + 2.0D0*C1/C2
      FA = CDABS(C3)
      RETURN

  380 T = Y(41+KA)
      BETA = Y(82+KA)
      X1 = (X(1)+ (1D0+X(2))*T)**2 + ((1D0-X(2))*BETA)**2
      X2 = (X(3)+ (1D0+X(4))*T)**2 + ((1D0-X(4))*BETA)**2
      X3 = (X(5)+ (1D0+X(6))*T)**2 + ((1D0-X(6))*BETA)**2
      X4 = (X(7)+ (1D0+X(8))*T)**2 + ((1D0-X(8))*BETA)**2
      IF (X2.EQ.0D0) X2 = 1D-30
      IF (X4.EQ.0D0) X4 = 1D-30
      FA = X(9)*SQRT(X1/X2)*SQRT(X3/X4) - ABS(1D0-2D0*Y(KA))
      RETURN

  390 FA = (X(1)-1.0D1)**2 + 5.0D0* (X(2)-1.2D1)**2 + X(3)**4 +
     +     3.0D0* (X(4)-1.1D1)**2 + 1.0D1*X(5)**6 + 7.0D0*X(6)**2 +
     +     X(7)**4 - 4.0D0*X(6)*X(7) - 1.0D1*X(6) - 8.0D0*X(7)
  400 GO TO (320,410,420,430,440) KA

  410 FA = FA + 1.0D1* (2.0D0*X(1)**2+3.0D0*X(2)**4+X(3)+4.0D0*X(4)**2+
     +     5.0D0*X(5)-1.27D2)
      RETURN

  420 FA = FA + 1.0D1* (7.0D0*X(1)+3.0D0*X(2)+1.0D1*X(3)**2+X(4)-X(5)-
     +     2.82D2)
      RETURN

  430 FA = FA + 1.0D1* (2.3D1*X(1)+X(2)**2+6.0D0*X(6)**2-8.0D0*X(7)-
     +     1.96D2)
      RETURN

  440 FA = FA + 1.0D1* (4.0D0*X(1)**2+X(2)**2-3.0D0*X(1)*X(2)+
     +     2.0D0*X(3)**2+5.0D0*X(6)-1.1D1*X(7))
      RETURN

  450 FA = X(1)**2 + X(2)**2 + X(1)*X(2) - 1.4D1*X(1) - 1.6D1*X(2) +
     +     (X(3)-1.0D1)**2 + 4.0D0* (X(4)-5.0D0)**2 + (X(5)-3.0D0)**2 +
     +     2.0D0* (X(6)-1.0D0)**2 + 5.0D0*X(7)**2 +
     +     7.0D0* (X(8)-1.1D1)**2 + 2.0D0* (X(9)-1.0D1)**2 +
     +     (X(10)-7.0D0)**2 + 4.5D1
  460 GO TO (320,470,480,490,500,510,520,530,540) KA

  470 FA = FA + 1.0D1* (3.0D0* (X(1)-2.0D0)**2+4.0D0* (X(2)-3.0D0)**2+
     +     2.0D0*X(3)**2-7.0D0*X(4)-1.2D2)
      RETURN

  480 FA = FA + 1.0D1* (5.0D0*X(1)**2+8.0D0*X(2)+ (X(3)-6.0D0)**2-
     +     2.0D0*X(4)-4.0D1)
      RETURN

  490 FA = FA + 1.0D1* (0.5D0* (X(1)-8.0D0)**2+2.0D0* (X(2)-4.0D0)**2+
     +     3.0D0*X(5)**2-X(6)-3.0D1)
      RETURN

  500 FA = FA + 1.0D1* (X(1)**2+2.0D0* (X(2)-2.0D0)**2-2.0D0*X(1)*X(2)+
     +     1.4D1*X(5)-6.0D0*X(6))
      RETURN

  510 FA = FA + 1.0D1* (4.0D0*X(1)+5.0D0*X(2)-3.0D0*X(7)+9.0D0*X(8)-
     +     1.05D2)
      RETURN

  520 FA = FA + 1.0D1* (1.0D1*X(1)-8.0D0*X(2)-1.7D1*X(7)+2.0D0*X(8))
      RETURN

  530 FA = FA + 1.0D1* (6.0D0*X(2)-3.0D0*X(1)+1.2D1* (X(9)-8.0D0)**2-
     +     7.0D0*X(10))
      RETURN

  540 FA = FA + 1.0D1* (2.0D0*X(2)-8.0D0*X(1)+5.0D0*X(9)-2.0D0*X(10)-
     +     1.2D1)
      RETURN

  550 FA = X(1)**2 + X(2)**2 + X(1)*X(2) - 1.4D1*X(1) - 1.6D1*X(2) +
     +     (X(3)-1.0D1)**2 + 4.0D0* (X(4)-5.0D0)**2 + (X(5)-3.0D0)**2 +
     +     2.0D0* (X(6)-1.0D0)**2 + 5.0D0*X(7)**2 +
     +     7.0D0* (X(8)-1.1D1)**2 + 2.0D0* (X(9)-1.0D1)**2 +
     +     (X(10)-7.0D0)**2 + (X(11)-9.0D0)**2 +
     +     1.0D1* (X(12)-1.0D0)**2 + 5.0D0* (X(13)-7.0D0)**2 +
     +     4.0D0* (X(14)-1.4D1)**2 + 2.7D1* (X(15)-1.0D0)**2 +
     +     X(16)**4 + (X(17)-2.0D0)**2 + 1.3D1* (X(18)-2.0D0)**2 +
     +     (X(19)-3.D0)**2 + X(20)**2 + 9.5D1
  560 GO TO (320,470,480,490,500,510,520,530,540,570,580,590,
     +       600,610,620,630,640,650) KA

  570 FA = FA + 1.0D1* (X(1)+X(2)+4.0D0*X(11)-2.1D1*X(12))
      RETURN

  580 FA = FA + 1.0D1* (X(1)**2+1.5D1*X(11)-8.0D0*X(12)-2.8D1)
      RETURN

  590 FA = FA + 1.0D1* (4.0D0*X(1)+9.0D0*X(2)+5.0D0*X(13)**2-
     +     9.0D0*X(14)-8.7D1)
      RETURN

  600 FA = FA + 1.0D1* (3.0D0*X(1)+4.0D0*X(2)+3.0D0* (X(13)-6.0D0)**2-
     +     1.4D1*X(14)-1.0D1)
      RETURN

  610 FA = FA + 1.0D1* (1.4D1*X(1)**2+3.5D1*X(15)-7.9D1*X(16)-9.2D1)
      RETURN

  620 FA = FA + 1.0D1* (1.5D1*X(2)**2+1.1D1*X(15)-6.1D1*X(16)-5.4D1)
      RETURN

  630 FA = FA + 1.0D1* (5.0D0*X(1)**2+2.0D0*X(2)+9.0D0*X(17)**4-X(18)-
     +     6.8D1)
      RETURN

  640 FA = FA + 1.0D1* (X(1)**2-X(2)+1.9D1*X(19)-2.0D1*X(20)+1.9D1)
      RETURN

  650 FA = FA + 1.0D1* (7.0D0*X(1)**2+5.0D0*X(2)**2+X(19)**2-
     +     3.0D1*X(20))
      RETURN

  660 X1 = 0.0D0
      DO 670 I = 1,N
          X3 = 1.0D0
          X4 = X(I)
          IF (I.EQ.1) X3 = 1.0D-8
          IF (I.EQ.4) X3 = 4.0D0
          IF (I.EQ.2 .AND. KA.EQ.1) X4 = X(I) + 2.0D0
          IF (I.EQ.2 .AND. KA.EQ.2) X4 = X(I) - 2.0D0
          X1 = X1 + X3*X4**2
  670 CONTINUE
      FA = EXP(X1)
      RETURN

  680 FA = 0.0D0
      DO 690 I = 1,N
          X1 = 1.0D0*DBLE(I+KA-1)
          X2 = X(I) - SIN(DBLE(2*I+KA-3))
          FA = FA + X1*EXP(X2**2)
  690 CONTINUE
      RETURN

  700 IF (KA.EQ.1) THEN
          FA = X(1)

      ELSE IF (KA.EQ.2) THEN
          FA = X(2) - X(1)**2 - 1.0D0

      ELSE
          T = DBLE(KA-2)/2.9D1
          X1 = 0.0D0
          X2 = X(1)
          DO 710 I = 2,N
              X1 = X1 + DBLE(I-1)*X(I)*T** (I-2)
              X2 = X2 + X(I)*T** (I-1)
  710     CONTINUE
          FA = X1 - X2**2 - 1.0D0
      END IF

      RETURN

  720 T = 1.0D-1*DBLE(KA-1)
      FA = Y(KA) - X(1)*EXP(-X(5)*T) - X(2)*EXP(-X(6)* (T-X(9))**2) -
     +     X(3)*EXP(-X(7)* (T-X(10))**2) - X(4)*EXP(-X(8)* (T-X(11))**2)
      RETURN

      END
* SUBROUTINE TAGU06             ALL SYSTEMS                 99/12/01
C PORTABILITY : ALL SYSTEMS
C 90/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
*  GRADIENTS OF PARTIAL FUNCTIONS IN THE MINIMAX CRITERION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  II  KA  INDEX OF THE PARTIAL FUNCTION.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  GA(N)  GRADIENT OF THE PARTIAL FUNCTION AT THE
*          SELECTED POINT.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*
      SUBROUTINE TAGU06(N,KA,X,GA,NEXT)
C     .. Parameters ..
      DOUBLE PRECISION PI
      PARAMETER (PI=3.14159265358979323846D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER KA,N,NEXT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION GA(N),X(N)
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION Y(123)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX C1,C2,C3
      DOUBLE PRECISION BETA,FA,T,X1,X2,X3,X4,X5,X6,X7,X8
      INTEGER I,J
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX CA(4),CB(4),CC(6)
      DOUBLE PRECISION XA(3),XB(3),XC(3)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,CDABS,CMPLX,CONJG,COS,DBLE,EXP,MIN,SIN,SQRT
C     ..
C     .. Common blocks ..
      COMMON /EMPR06/Y
C     ..
      GO TO (10,50,60,70,140,190,240,250,260,
     +       270,280,290,300,310,330,340,350,400,
     +       410,470,570,680,710,730,770) NEXT

   10 X1 = X(1)*X(1)
      X2 = X(2)*X(2)
      X3 = X(1) + X(1)
      X4 = X(2) + X(2)
      GO TO (20,30,40) KA

   20 GA(1) = X3
      GA(2) = (X2+X2)*X4
      RETURN

   30 GA(1) = -4.0D0 + X3
      GA(2) = -4.0D0 + X4
      RETURN

   40 FA = 2.0D0*EXP(X(2)-X(1))
      GA(1) = -FA
      GA(2) = +FA
      RETURN

   50 X1 = 1.0D0/ (X(1)+1.0D-1)**2
      GA(2) = 2.0D0*X(2)
      IF (KA.EQ.1) THEN
          GA(1) = 0.5D0* (1.0D0+X1)

      ELSE IF (KA.EQ.2) THEN
          GA(1) = 0.5D0* (-1.0D0+X1)

      ELSE IF (KA.EQ.3) THEN
          GA(1) = 0.5D0* (1.0D0-X1)
      END IF

      RETURN

   60 X1 = X(1)**2 + X(2)**2
      X2 = SQRT(X1)
      X3 = COS(X2)
      X4 = SIN(X2)
      IF (KA.EQ.1) THEN
          X5 = 2.0D0* (X(1)-X2*X3)
          X6 = - (X3/X2-X4)
          GA(1) = X5* (X(1)*X6+1.0D0) + 1.0D-2*X(1)
          GA(2) = X5*X(2)*X6 + 1.0D-2*X(2)

      ELSE IF (KA.EQ.2) THEN
          X5 = 2.0D0* (X(2)-X2*X4)
          X6 = - (X4/X2+X3)
          GA(1) = X5*X(1)*X6 + 1.0D-2*X(1)
          GA(2) = X5* (X(2)*X6+1.0D0) + 1.0D-2*X(2)
      END IF

      RETURN

   70 GO TO (80,90,100,110,120,130) KA

   80 GA(1) = 2.0D0*X(1)
      GA(2) = 2.0D0*X(2)
      GA(3) = 2.0D0*X(3)
      RETURN

   90 GA(1) = 2.0D0*X(1)
      GA(2) = 2.0D0*X(2)
      GA(3) = 2.0D0* (X(3)-2.0D0)
      RETURN

  100 GA(1) = 1.0D0
      GA(2) = 1.0D0
      GA(3) = 1.0D0
      RETURN

  110 GA(1) = 1.0D0
      GA(2) = 1.0D0
      GA(3) = -1.0D0
      RETURN

  120 GA(1) = 6.0D0*X(1)**2 - 4.0D0* (5.0D0*X(3)-X(1)+1.0D0)
      GA(2) = 1.2D1*X(2)
      GA(3) = 2.0D1* (5.0D0*X(3)-X(1)+1.0D0)
      RETURN

  130 GA(1) = 2.0D0*X(1)
      GA(2) = 0.0D0
      GA(3) = -9.0D0
      RETURN

  140 X1 = X(1)*X(1)
      X2 = X(2)*X(2)
      X3 = X(3)*X(3)
      X4 = X(4)*X(4)
      X5 = X(1) + X(1)
      X6 = X(2) + X(2)
      X7 = X(3) + X(3)
      X8 = X(4) + X(4)
      GA(1) = X5 - 5.0D0
      GA(2) = X6 - 5.0D0
      GA(3) = X7 + X7 - 2.1D1
      GA(4) = X8 + 7.0D0
  150 GO TO (320,160,170,180) KA

  160 GA(1) = GA(1) + 1.0D1* (X5+1.0D0)
      GA(2) = GA(2) + 1.0D1* (X6-1.0D0)
      GA(3) = GA(3) + 1.0D1* (X7+1.0D0)
      GA(4) = GA(4) + 1.0D1* (X8-1.0D0)
      RETURN

  170 GA(1) = GA(1) + 1.0D1* (X5-1.0D0)
      GA(2) = GA(2) + 1.0D1* (X6+X6)
      GA(3) = GA(3) + 1.0D1*X7
      GA(4) = GA(4) + 1.0D1* (X8+X8-1.0D0)
      RETURN

  180 GA(1) = GA(1) + 1.0D1* (X5+2.0D0)
      GA(2) = GA(2) + 1.0D1* (X6-1.0D0)
      GA(3) = GA(3) + 1.0D1*X7
      GA(4) = GA(4) - 1.0D1
      RETURN

  190 X1 = X(1) - (X(4)+1.0D0)**4
      X2 = X1*X1
      X3 = X(2) - X2*X2
      X4 = X1*X3
      X5 = -4.0D0* (X(4)+1.0D0)**3
      GA(1) = 2.0D0*X1 - 8.0D0*X4 - 5.0D0* (1.0D0-4.0D0*X1)
      GA(2) = 2.0D0*X3 - 5.0D0
      GA(3) = 4.0D0*X(3) - 2.1D1
      GA(4) = 2.0D0*X1*X5 - 8.0D0*X4*X5 + 2.0D0*X(4) -
     +        5.0D0* (X5-4.0D0*X1*X5) + 7.0D0
      GO TO (200,210,220,230) KA

  200 CONTINUE
      RETURN

  210 GA(1) = GA(1) + 1.0D1* (2.0D0*X1-8.0D0*X4+1.0D0+4.0D0*X1)
      GA(2) = GA(2) + 1.0D1* (2.0D0*X3-1.0D0)
      GA(3) = GA(3) + 1.0D1* (2.0D0*X(3)+1.0D0)
      GA(4) = GA(4) + 1.0D1* (2.0D0*X1*X5-8.0D0*X4*X5+2.0D0*X(4)+X5+
     +        4.0D0*X1*X5-1.0D0)
      RETURN

  220 GA(1) = GA(1) + 1.0D1* (2.0D0*X1-1.6D1*X4-1.0D0)
      GA(2) = GA(2) + 1.0D1* (4.0D0*X3)
      GA(3) = GA(3) + 1.0D1* (2.0D0*X(3))
      GA(4) = GA(4) + 1.0D1* (2.0D0*X1*X5-1.6D1*X4*X5+4.0D0*X(4)-X5-
     +        1.0D0)
      RETURN

  230 GA(1) = GA(1) + 1.0D1* (2.0D0*X1-8.0D0*X4+2.0D0+4.0D0*X1)
      GA(2) = GA(2) + 1.0D1* (2.0D0*X3-1.0D0)
      GA(3) = GA(3) + 1.0D1* (2.0D0*X(3))
      GA(4) = GA(4) + 1.0D1* (2.0D0*X1*X5-8.0D0*X4*X5+2.0D0*X5+
     +        4.0D0*X1*X5-1.0D0)
      RETURN

  240 T = Y(KA)
      X1 = EXP(-X(1)*T)/X(2)
      X2 = SIN(X(2)*T)
      X3 = COS(X(2)*T)
      GA(1) = -T*X(3)*X1*X2
      GA(2) = X(3)*X1* (T*X3-X2/X(2))
      GA(3) = X1*X2
      RETURN

  250 C1 = DBLE(16-KA)
      C2 = DBLE(MIN(KA,16-KA))
      C3 = DBLE(KA)/ (C1*X(2)+C2*X(3))**2
      GA(1) = -1.0D0
      GA(2) = C1*C3
      GA(3) = C2*C3
      RETURN

  260 T = Y(KA+11)
      X1 = X(1)*T* (T+X(2))
      X2 = ((T+X(3))*T+X(4))
      X3 = X1/X2**2
      GA(1) = -T* (T+X(2))/X2
      GA(2) = -T*X(1)/X2
      GA(3) = T*X3
      GA(4) = X3
      RETURN

  270 T = 0.2D0*DBLE(KA)
      GA(1) = 2.0D0* (X(1)+X(2)*T-EXP(T))
      GA(2) = 2.0D0* (X(1)+X(2)*T-EXP(T))*T
      GA(3) = 2.0D0* (X(3)+X(4)*SIN(T)-COS(T))
      GA(4) = 2.0D0* (X(3)+X(4)*SIN(T)-COS(T))*SIN(T)
      RETURN

  280 T = Y(KA)
      X1 = -2.0D0* ((X(1)*T+X(2))*T+X(3))
      GA(1) = X1*T**2
      GA(2) = X1*T
      GA(3) = X1
      GA(4) = 1.0D0
      RETURN

  290 T = Y(KA)
      X1 = EXP(X(3)*T)
      X2 = EXP(X(4)*T)
      GA(1) = X1
      GA(2) = X2
      GA(3) = X(1)*T*X1
      GA(4) = X(2)*T*X2
      RETURN

  300 T = Y(KA)
      X1 = T + X(2) + 1.0D0/ (X(3)*T+X(4))
      IF (X1.EQ.0D0) X1 = 1D-30
      X2 = X1/ ((T+1.0D0)*Y(61+KA))
      GA(1) = ABS(X2)** (T+5.0D-1)
      GA(2) = X(1)*GA(1)* (T+5.0D-1)/X1
      GA(4) = -GA(2)/ (X(3)*T+X(4))**2
      GA(3) = GA(4)*T
      RETURN

  310 T = 0.1D0*DBLE(KA-1) - 1.0D0
      X1 = X(1) + T*X(2)
      X2 = 1.0D0/ (1.0D0+T* (X(3)+T* (X(4)+T*X(5))))
      X3 = X1*X2 - EXP(T)
      GA(1) = X2
      GA(2) = X2*T
      GA(3) = -X1*X2*X2*T
      GA(4) = GA(3)*T
      GA(5) = GA(4)*T
  320 RETURN

  330 T = Y(KA)
      X1 = 1.0D0/ (1.0D0+T* (X(4)+T*X(5)))
      X2 = X(1) + T* (X(2)+T*X(3))
      GA(1) = X1
      GA(2) = X1*T
      GA(3) = X1*T*T
      GA(4) = -X2*X1**2*T
      GA(5) = -X2*X1**2*T*T
      RETURN

  340 T = 0.1D0*DBLE(KA-1)
      X1 = EXP(-X(2)*T)
      X2 = COS(X(3)*T+X(4))
      X3 = SIN(X(3)*T+X(4))
      X4 = EXP(-X(6)*T)
      GA(1) = X1*X2
      GA(2) = -X1*X2*X(1)*T
      GA(3) = -X1*X3*X(1)*T
      GA(4) = -X1*X3*X(1)
      GA(5) = X4
      GA(6) = -X4*X(5)*T
      RETURN

  350 BETA = 0.5D0*PI*Y(KA)
      DO 360 I = 1,3
          J = I + I
          XA(I) = X(J-1)
          XB(I) = X(J)
  360 CONTINUE
      CA(4) = CMPLX(1.0D0,0.0D0)
      CB(4) = 1.0D1*CA(4)
      DO 370 J = 1,3
          I = 4 - J
          XC(I) = BETA*XA(I)
          T = XC(I)
          X1 = COS(T)
          X2 = SIN(T)
          C1 = CMPLX(X1,0.0D0)
          C2 = CMPLX(0.0D0, (X2*XB(I)))
          C3 = CMPLX(0.0D0, (X2/XB(I)))
          CB(I) = C1*CB(I+1) + C2*CA(I+1)
          CA(I) = C3*CB(I+1) + C1*CA(I+1)
  370 CONTINUE
      C1 = -CA(1)
      C2 = CB(1) - C1
      C3 = 1.0D0 + 2.0D0*C1/C2
      FA = CDABS(C3)
      C3 = CONJG(C3)
      C1 = 2.0D0/C2
      DO 380 I = 1,3
          T = XC(I)
          J = I + I
          CC(J) = (CB(I)*CA(I)-CB(I+1)*CA(I+1))/ (C2*XB(I))
          CC(J-1) = BETA* (CB(I)*CA(I+1)-CB(I+1)*CA(I))/ (C2*SIN(T))
  380 CONTINUE
      DO 390 I = 1,6
          GA(I) = DBLE(C1*C3*CC(I))/FA
  390 CONTINUE
      RETURN

  400 T = Y(41+KA)
      BETA = Y(82+KA)
      X1 = (X(1)+ (1D0+X(2))*T)**2 + ((1D0-X(2))*BETA)**2
      X2 = (X(3)+ (1D0+X(4))*T)**2 + ((1D0-X(4))*BETA)**2
      X3 = (X(5)+ (1D0+X(6))*T)**2 + ((1D0-X(6))*BETA)**2
      X4 = (X(7)+ (1D0+X(8))*T)**2 + ((1D0-X(8))*BETA)**2
      IF (X1.EQ.0D0) X1 = 1D-30
      IF (X2.EQ.0D0) X2 = 1D-30
      IF (X3.EQ.0D0) X3 = 1D-30
      IF (X4.EQ.0D0) X4 = 1D-30
      FA = SQRT(X1/X2)*SQRT(X3/X4)
      GA(9) = FA
      FA = X(9)*FA
      GA(1) = FA/X1* (X(1)+T* (1D0+X(2)))
      GA(2) = FA/X1* (X(2)+2D0*T*T-1D0+X(1)*T)
      GA(3) = -FA/X2* (X(3)+T* (1D0+X(4)))
      GA(4) = -FA/X2* (X(4)+2D0*T*T-1D0+X(3)*T)
      GA(5) = FA/X3* (X(5)+T* (1D0+X(6)))
      GA(6) = FA/X3* (X(6)+2D0*T*T-1D0+X(5)*T)
      GA(7) = -FA/X4* (X(7)+T* (1D0+X(8)))
      GA(8) = -FA/X4* (X(8)+2D0*T*T-1D0+X(7)*T)
      RETURN

  410 GA(1) = 2.0D0* (X(1)-1.0D1)
      GA(2) = 1.0D1* (X(2)-1.2D1)
      GA(3) = 4.0D0*X(3)**3
      GA(4) = 6.0D0* (X(4)-1.1D1)
      GA(5) = 6.0D1*X(5)**5
      GA(6) = 1.4D1*X(6) - 4.0D0*X(7) - 1.0D1
      GA(7) = 4.0D0*X(7)**3 - 4.0D0*X(6) - 8.0D0
  420 GO TO (320,430,440,450,460) KA

  430 GA(1) = GA(1) + 4.0D1*X(1)
      GA(2) = GA(2) + 1.2D2*X(2)**3
      GA(3) = GA(3) + 1.0D1
      GA(4) = GA(4) + 8.0D1*X(4)
      GA(5) = GA(5) + 5.0D1
      RETURN

  440 GA(1) = GA(1) + 7.0D1
      GA(2) = GA(2) + 3.0D1
      GA(3) = GA(3) + 2.0D2*X(3)
      GA(4) = GA(4) + 1.0D1
      GA(5) = GA(5) - 1.0D1
      RETURN

  450 GA(1) = GA(1) + 2.3D2
      GA(2) = GA(2) + 2.0D1*X(2)
      GA(6) = GA(6) + 1.2D2*X(6)
      GA(7) = GA(7) - 8.0D1
      RETURN

  460 GA(1) = GA(1) + 8.0D1*X(1) - 3.0D1*X(2)
      GA(2) = GA(2) + 2.0D1*X(2) - 3.0D1*X(1)
      GA(3) = GA(3) + 4.0D1*X(3)
      GA(6) = GA(6) + 5.0D1
      GA(7) = GA(7) - 1.1D2
      RETURN

  470 GA(1) = 2.0D0*X(1) + X(2) - 1.4D1
      GA(2) = 2.0D0*X(2) + X(1) - 1.6D1
      GA(3) = 2.0D0* (X(3)-1.0D1)
      GA(4) = 8.0D0* (X(4)-5.0D0)
      GA(5) = 2.0D0* (X(5)-3.0D0)
      GA(6) = 4.0D0* (X(6)-1.0D0)
      GA(7) = 1.0D1*X(7)
      GA(8) = 1.4D1* (X(8)-1.1D1)
      GA(9) = 4.0D0* (X(9)-1.0D1)
      GA(10) = 2.0D0* (X(10)-7.0D0)
  480 GO TO (320,490,500,510,520,530,540,550,560) KA

  490 GA(1) = GA(1) + 6.0D1* (X(1)-2.0D0)
      GA(2) = GA(2) + 8.0D1* (X(2)-3.0D0)
      GA(3) = GA(3) + 4.0D1*X(3)
      GA(4) = GA(4) - 7.0D1
      RETURN

  500 GA(1) = GA(1) + 1.0D2*X(1)
      GA(2) = GA(2) + 8.0D1
      GA(3) = GA(3) + 2.0D1* (X(3)-6.0D0)
      GA(4) = GA(4) - 2.0D1
      RETURN

  510 GA(1) = GA(1) + 1.0D1* (X(1)-8.0D0)
      GA(2) = GA(2) + 4.0D1* (X(2)-4.0D0)
      GA(5) = GA(5) + 6.0D1*X(5)
      GA(6) = GA(6) - 1.0D1
      RETURN

  520 GA(1) = GA(1) + 2.0D1*X(1) - 2.0D1*X(2)
      GA(2) = GA(2) + 4.0D1* (X(2)-2.0D0) - 2.0D1*X(1)
      GA(5) = GA(5) + 1.4D2
      GA(6) = GA(6) - 6.0D1
      RETURN

  530 GA(1) = GA(1) + 4.0D1
      GA(2) = GA(2) + 5.0D1
      GA(7) = GA(7) - 3.0D1
      GA(8) = GA(8) + 9.0D1
      RETURN

  540 GA(1) = GA(1) + 1.0D2
      GA(2) = GA(2) - 8.0D1
      GA(7) = GA(7) - 1.7D2
      GA(8) = GA(8) + 2.0D1
      RETURN

  550 GA(1) = GA(1) - 3.0D1
      GA(2) = GA(2) + 6.0D1
      GA(9) = GA(9) + 2.4D2* (X(9)-8.0D0)
      GA(10) = GA(10) - 7.0D1
      RETURN

  560 GA(1) = GA(1) - 8.0D1
      GA(2) = GA(2) + 2.0D1
      GA(9) = GA(9) + 5.0D1
      GA(10) = GA(10) - 2.0D1
      RETURN

  570 GA(1) = 2.0D0*X(1) + X(2) - 1.4D1
      GA(2) = 2.0D0*X(2) + X(1) - 1.6D1
      GA(3) = 2.0D0* (X(3)-1.0D1)
      GA(4) = 8.0D0* (X(4)-5.0D0)
      GA(5) = 2.0D0* (X(5)-3.0D0)
      GA(6) = 4.0D0* (X(6)-1.0D0)
      GA(7) = 1.0D1*X(7)
      GA(8) = 1.4D1* (X(8)-1.1D1)
      GA(9) = 4.0D0* (X(9)-1.0D1)
      GA(10) = 2.0D0* (X(10)-7.0D0)
      GA(11) = 2.0D0* (X(11)-9.0D0)
      GA(12) = 2.0D1* (X(12)-1.0D0)
      GA(13) = 1.0D1* (X(13)-7.0D0)
      GA(14) = 8.0D0* (X(14)-1.4D1)
      GA(15) = 5.4D1* (X(15)-1.0D0)
      GA(16) = 4.0D0*X(16)**3
      GA(17) = 2.0D0* (X(17)-2.0D0)
      GA(18) = 2.6D1* (X(18)-2.0D0)
      GA(19) = 2.0D0* (X(19)-3.0D0)
      GA(20) = 2.0D0*X(20)
  580 GO TO (320,490,500,510,520,530,540,550,560,590,600,610,
     +       620,630,640,650,660,670) KA

  590 GA(1) = GA(1) + 1.0D1
      GA(2) = GA(2) + 1.0D1
      GA(11) = GA(11) + 4.0D1
      GA(12) = GA(12) - 2.1D2
      RETURN

  600 GA(1) = GA(1) + 2.0D1*X(1)
      GA(11) = GA(11) + 1.5D2
      GA(12) = GA(12) - 8.0D1
      RETURN

  610 GA(1) = GA(1) + 4.0D1
      GA(2) = GA(2) + 9.0D1
      GA(13) = GA(13) + 1.0D2*X(13)
      GA(14) = GA(14) - 9.0D1
      RETURN

  620 GA(1) = GA(1) + 3.0D1
      GA(2) = GA(2) + 4.0D1
      GA(13) = GA(13) + 6.0D1* (X(13)-6.0D0)
      GA(14) = GA(14) - 1.4D2
      RETURN

  630 GA(1) = GA(1) + 2.8D2*X(1)
      GA(15) = GA(15) + 3.5D2
      GA(16) = GA(16) - 7.9D2
      RETURN

  640 GA(2) = GA(2) + 3.0D2*X(2)
      GA(15) = GA(15) + 1.1D2
      GA(16) = GA(16) - 6.1D2
      RETURN

  650 GA(1) = GA(1) + 1.0D2*X(1)
      GA(2) = GA(2) + 2.0D1
      GA(17) = GA(17) + 3.6D2*X(17)**3
      GA(18) = GA(18) - 1.0D1
      RETURN

  660 GA(1) = GA(1) + 2.0D1*X(1)
      GA(2) = GA(2) - 1.0D1
      GA(19) = GA(19) + 1.9D2
      GA(20) = GA(20) - 2.0D2
      RETURN

  670 GA(1) = GA(1) + 1.4D2*X(1)
      GA(2) = GA(2) + 1.0D2*X(2)
      GA(19) = GA(19) + 2.0D1*X(19)
      GA(20) = GA(20) - 3.0D2
      RETURN

  680 X1 = 0.0D0
      DO 690 I = 1,N
          X3 = 1.0D0
          X4 = X(I)
          IF (I.EQ.1) X3 = 1.0D-8
          IF (I.EQ.4) X3 = 4.0D0
          IF (I.EQ.2 .AND. KA.EQ.1) X4 = X(I) + 2.0D0
          IF (I.EQ.2 .AND. KA.EQ.2) X4 = X(I) - 2.0D0
          X1 = X1 + X3*X4**2
  690 CONTINUE
      X2 = EXP(X1)
      DO 700 I = 1,N
          X3 = 2.0D0
          X4 = X(I)
          IF (I.EQ.1) X3 = 2.0D-8
          IF (I.EQ.4) X3 = 8.0D0
          IF (I.EQ.2 .AND. KA.EQ.1) X4 = X(I) + 2.0D0
          IF (I.EQ.2 .AND. KA.EQ.2) X4 = X(I) - 2.0D0
          GA(I) = X2*X3*X4
  700 CONTINUE
      RETURN

  710 DO 720 I = 1,N
          X1 = 1.0D0*DBLE(I+KA-1)
          X2 = X(I) - SIN(DBLE(2*I+KA-3))
          GA(I) = 2.0D0*X1*X2*EXP(X2**2)
  720 CONTINUE
      RETURN

  730 IF (KA.LE.2) THEN
          DO 740 I = 2,N
              GA(I) = 0.0D0
  740     CONTINUE
          IF (KA.EQ.1) THEN
              GA(1) = 1.0D0

          ELSE IF (KA.EQ.2) THEN
              GA(1) = -2.0D0*X(1)
              GA(2) = 1.0D0
          END IF

      ELSE
          GA(1) = 0.0D0
          T = DBLE(KA-2)/2.9D1
          X2 = X(1)
          DO 750 I = 2,N
              X2 = X2 + X(I)*T** (I-1)
  750     CONTINUE
          DO 760 I = 1,N
              IF (I.GT.1) GA(I) = DBLE(I-1)*T** (I-2)
              GA(I) = GA(I) - 2.0D0*X2*T** (I-1)
  760     CONTINUE
      END IF

      RETURN

  770 T = 1.0D-1*DBLE(KA-1)
      X1 = EXP(-X(5)*T)
      X2 = EXP(-X(6)* (T-X(9))**2)
      X3 = EXP(-X(7)* (T-X(10))**2)
      X4 = EXP(-X(8)* (T-X(11))**2)
      GA(1) = -X1
      GA(2) = -X2
      GA(3) = -X3
      GA(4) = -X4
      GA(5) = X1*X(1)*T
      GA(6) = X2*X(2)* (T-X(9))**2
      GA(7) = X3*X(3)* (T-X(10))**2
      GA(8) = X4*X(4)* (T-X(11))**2
      GA(9) = -2.0D0*X2*X(2)*X(6)* (T-X(9))
      GA(10) = -2.0D0*X3*X(3)*X(7)* (T-X(10))
      GA(11) = -2.0D0*X4*X(4)*X(8)* (T-X(11))
      RETURN

      END
* SUBROUTINE TAHD06             ALL SYSTEMS                 99/12/01
C PORTABILITY : ALL SYSTEMS
C 95/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
*  HESSIAN MATRICES OF PARTIAL FUNCTIONS IN THE MINIMAX CRITERION.
*  DENSE VERSION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  II  KA  INDEX OF THE PARTIAL FUNCTION.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  HA(N*(N+1)/2)  HESSIAN MATRIX OF THE PARTIAL FUNCTION
*         AT THE SELECTED POINT.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*
      SUBROUTINE TAHD06(N,KA,X,HA,NEXT)
C     .. Parameters ..
      DOUBLE PRECISION PI
      PARAMETER (PI=3.14159265358979323846D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER KA,N,NEXT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION HA(N* (N+1)/2),X(N)
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION Y(123)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX C1,C2,C3,CI,S1,S2,S3,S4
      DOUBLE PRECISION BETA,FA,T,X1,X2,X3,X4,X5,X6,X7,X8
      INTEGER I,J,L
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX CA(4),CB(4),CC(6),DD(6)
      DOUBLE PRECISION CT(3),GA(8),ST(3),XA(3),XB(3),XC(3)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,CDABS,CMPLX,CONJG,COS,DBLE,EXP,MIN,SIN,SQRT
C     ..
C     .. Common blocks ..
      COMMON /EMPR06/Y
C     ..
      GO TO (10,20,30,40,100,150,200,210,230,
     +       240,250,270,290,300,310,330,350,430,
     +       470,540,620,720,760,790,830) NEXT

   10 HA(1) = 2.0D0
      HA(2) = 0.0D0
      HA(3) = 2.0D0
      IF (KA.EQ.1) HA(3) = 12.0D0*X(2)*X(2)
      IF (KA.EQ.3) THEN
          HA(1) = 2.0D0*EXP(X(2)-X(1))
          HA(2) = -HA(1)
          HA(3) = HA(1)
      END IF

      RETURN

   20 CONTINUE
      HA(2) = 0.0D0
      HA(3) = 2.0D0
      X1 = 1.0D0/ (X(1)+1.0D-1)**3
      IF (KA.EQ.3) THEN
          HA(1) = X1

      ELSE
          HA(1) = -X1
      END IF

      RETURN

   30 X1 = X(1)**2 + X(2)**2
      X2 = SQRT(X1)
      X3 = COS(X2)
      X4 = SIN(X2)
      IF (KA.EQ.1) THEN
          X5 = X(1) - X2*X3
          X6 = - (X3/X2-X4)
          X7 = X2 + 1.0D0/X2
          X7 = (X7*X3+X4)/X2**2
          X1 = 1.0D0 + X6*X(1)
          X2 = X6*X(2)
          HA(1) = 2.0D0* (X1**2+X5* (X6+X7*X(1)**2)) + 1.0D-2
          HA(2) = 2.0D0*X(2)* (X6*X1+X5*X7*X(1))
          HA(3) = 2.0D0* (X2**2+X5* (X6+X7*X(2)**2)) + 1.0D-2

      ELSE IF (KA.EQ.2) THEN
          X5 = X(2) - X2*X4
          X6 = - (X4/X2+X3)
          X7 = X2 + 1.0D0/X2
          X7 = (X7*X4-X3)/X2**2
          X1 = X6*X(1)
          X2 = 1.0D0 + X6*X(2)
          HA(1) = 2.0D0* (X1**2+X5* (X6+X7*X(1)**2)) + 1.0D-2
          HA(2) = 2.0D0*X(1)* (X6*X2+X5*X7*X(2))
          HA(3) = 2.0D0* (X2**2+X5* (X6+X7*X(2)**2)) + 1.0D-2
      END IF

      RETURN

   40 DO 50 I = 1,N* (N+1)/2
          HA(I) = 0.0D0
   50 CONTINUE
      GO TO (60,60,70,70,80,90) KA

   60 HA(1) = 2.0D0
      HA(3) = 2.0D0
      HA(6) = 2.0D0
      RETURN

   70 RETURN

   80 HA(1) = 1.6D0
      HA(3) = 1.2D1
      HA(4) = -2.0D1
      HA(6) = 1.0D2
      RETURN

   90 HA(1) = 2.0D0
      RETURN

  100 DO 110 I = 1,N* (N+1)/2
          HA(I) = 0.0D0
  110 CONTINUE
      HA(1) = 2.0D0
      HA(3) = 2.0D0
      HA(6) = 4.0D0
      HA(10) = 2.0D0
      IF (KA.GT.1) THEN
          HA(1) = HA(1) + 2.0D1
          HA(3) = HA(3) + 2.0D1
          HA(6) = HA(6) + 2.0D1
          HA(10) = HA(10) + 2.0D1
          GO TO (140,120,130) KA - 1

  120     HA(3) = HA(3) + 2.0D1
          HA(10) = HA(10) + 2.0D1
          GO TO 140

  130     HA(10) = HA(10) - 2.0D1
  140     CONTINUE
      END IF

      RETURN

  150 X1 = X(1) - (X(4)+1.0D0)**4
      X2 = X1*X1
      X3 = X(2) - X2*X2
      X4 = X3*X3
      X5 = -4.0D0* (X(4)+1.0D0)**3
      X6 = -1.2D1* (X(4)+1.0D0)**2
      X7 = X1*X6 + X5*X5
      X8 = X3 - 4.0D0*X1*X1
      HA(1) = 2.0D0 - 8.0D0*X8 + 2.0D1
      HA(2) = -8.0D0*X1
      HA(3) = 2.0D0
      HA(4) = 0.0D0
      HA(5) = 0.0D0
      HA(6) = 4.0D0
      HA(7) = 2.0D0*X5 - 8.0D0*X5*X8 + 2.0D1*X5
      HA(8) = -8.0D0*X1*X5
      HA(9) = 0.0D0
      HA(10) = 2.0D0*X7 - 8.0D0* (X5*X5*X8+X1*X3*X6) -
     +         5.0D0* (X6-4.0D0*X7) + 2.0D0
      GO TO (160,170,180,190) KA

  160 CONTINUE
      RETURN

  170 CONTINUE
      HA(1) = HA(1) + 1.0D1* (2.0D0-8.0D0*X8+4.0D0)
      HA(2) = HA(2) - 8.0D1*X1
      HA(3) = HA(3) + 2.0D1
      HA(6) = HA(6) + 2.0D1
      HA(7) = HA(7) + 1.0D1* (2.0D0*X5-8.0D0*X5*X8+4.0D0*X5)
      HA(8) = HA(8) - 8.0D1*X1*X5
      HA(10) = HA(10) + 1.0D1* (2.0D0*X7-8.0D0* (X5*X5*X8+X1*X3*X6)+X6+
     +         4.0D0*X7+2.0D0)
      RETURN

  180 CONTINUE
      HA(1) = HA(1) + 1.0D1* (2.0D0-1.6D1*X8)
      HA(2) = HA(2) - 1.6D2*X1
      HA(3) = HA(3) + 4.0D1
      HA(6) = HA(6) + 2.0D1
      HA(7) = HA(7) + 1.0D1* (2.0D0*X5-1.6D1*X5*X8)
      HA(8) = HA(8) - 1.6D2*X1*X5
      HA(10) = HA(10) + 1.0D1* (2.0D0*X7-1.6D1* (X5*X5*X8+X1*X3*X6)-X6+
     +         4.0D0)
      RETURN

  190 CONTINUE
      HA(1) = HA(1) + 1.0D1* (2.0D0-8.0D0*X8+4.0D0)
      HA(2) = HA(2) - 8.0D1*X1
      HA(3) = HA(3) + 2.0D1
      HA(6) = HA(6) + 2.0D1
      HA(7) = HA(7) + 1.0D1* (2.0D0*X5-8.0D0*X5*X8+4.0D0*X5)
      HA(8) = HA(8) - 8.0D1*X1*X5
      HA(10) = HA(10) + 1.0D1* (2.0D0*X7-8.0D0* (X5*X5*X8+X1*X3*X6)+
     +         2.0D0*X6+4.0D0*X7)
      RETURN

  200 T = Y(KA)
      X1 = EXP(-X(1)*T)/X(2)
      X2 = SIN(X(2)*T)
      X3 = COS(X(2)*T)
      X4 = T*X3 - X2/X(2)
      HA(1) = X(3)*X1*X2*T**2
      HA(2) = -X(3)*X1*X4*T
      HA(3) = -X(3)*X1* (T**2*X2+2.0D0*T*X3/X(2)-2.0D0*X2/X(2)**2)
      HA(4) = -X1*X2*T
      HA(5) = X1*X4*T
      HA(6) = 0.0D0
      RETURN

  210 DO 220 I = 1,N* (N+1)/2
          HA(I) = 0.0D0
  220 CONTINUE
      C1 = DBLE(16-KA)
      C2 = DBLE(MIN(KA,16-KA))
      C3 = -2.0D0*DBLE(KA)/ (C1*X(2)+C2*X(3))**3
      HA(3) = C1*C1*C3
      HA(5) = C1*C2*C3
      HA(6) = C2*C2*C3
      RETURN

  230 T = Y(KA+11)
      X1 = X(1)*T* (T+X(2))
      X2 = ((T+X(3))*T+X(4))
      X3 = 1.0D0/X2**2
      X4 = T* (T+X(2))*X3
      X5 = -2.0D0*X1*X3/X2
      HA(1) = 0.0D0
      HA(2) = -T/X2
      HA(3) = 0.0D0
      HA(7) = X4
      HA(4) = T*X4
      HA(8) = T*X(1)*X3
      HA(5) = T*HA(8)
      HA(10) = X5
      HA(9) = T*X5
      HA(6) = T*HA(9)
      RETURN

  240 T = 0.2D0*DBLE(KA)
      HA(1) = 2.0D0
      HA(2) = 2.0D0*T
      HA(3) = 2.0D0*T*T
      HA(4) = 0.0D0
      HA(5) = 0.0D0
      HA(6) = 2.0D0
      HA(7) = 0.0D0
      HA(8) = 0.0D0
      HA(9) = 2.0D0*SIN(T)
      HA(10) = 2.0D0*SIN(T)**2
      RETURN

  250 T = Y(KA)
      DO 260 I = 7,10
          HA(I) = 0.0D0
  260 CONTINUE
      HA(6) = -2.0D0
      HA(5) = HA(6)*T
      HA(4) = HA(5)*T
      HA(3) = HA(4)
      HA(2) = HA(3)*T
      HA(1) = HA(2)*T
      RETURN

  270 T = Y(KA)
      DO 280 I = 1,10
          HA(I) = 0.0D0
  280 CONTINUE
      X1 = EXP(X(3)*T)
      X2 = EXP(X(4)*T)
      HA(4) = X1*T
      HA(8) = X2*T
      HA(6) = X(1)*X1*T**2
      HA(10) = X(2)*X2*T**2
      RETURN

  290 T = Y(KA)
      X1 = T + X(2) + 1.0D0/ (X(3)*T+X(4))
      IF (X1.EQ.0D0) X1 = 1D-30
      X2 = X1/ ((T+1.0D0)*Y(61+KA))
      X3 = X(3)*T + X(4)
      IF (X3.EQ.0D0) X3 = 1D-30
      HA(1) = 0.0D0
      HA(2) = ABS(X2)** (T+5.0D-1)* (T+5.0D-1)/X1
      BETA = X(1)*HA(2)/X1
      HA(3) = BETA* (T-5.0D-1)
      HA(7) = -HA(2)/ (X3*X3)
      HA(4) = HA(7)*T
      HA(8) = -HA(3)/ (X3*X3)
      HA(5) = HA(8)*T
      HA(10) = BETA* (T+1.5D0+ (X3+X3)* (T+X(2)))/X3**4
      HA(6) = HA(10)*T*T
      HA(9) = HA(10)*T
      RETURN

  300 T = 0.1D0*DBLE(KA-1) - 1.0D0
      HA(1) = 0.0D0
      HA(2) = 0.0D0
      HA(3) = 0.0D0
      X1 = X(1) + T*X(2)
      X2 = 1.0D0/ (1.0D0+T* (X(3)+T* (X(4)+T*X(5))))
      X3 = -T*X2*X2
      X4 = -2.0D0*T*X1*X2*X3
      HA(4) = X3
      HA(5) = HA(4)*T
      HA(7) = HA(5)
      HA(8) = HA(7)*T
      HA(11) = HA(8)
      HA(12) = HA(11)*T
      HA(6) = X4
      HA(9) = HA(6)*T
      HA(10) = HA(9)*T
      HA(13) = HA(10)
      HA(14) = HA(13)*T
      HA(15) = HA(14)*T
      RETURN

  310 T = Y(KA)
      DO 320 I = 1,6
          HA(I) = 0.0D0
  320 CONTINUE
      X1 = 1.0D0/ (1.0D0-T* (X(4)-T*X(5)))
      X2 = X(1) + T* (X(2)+T*X(3))
      X3 = X1*X1
      HA(7) = -T*X3
      HA(8) = HA(7)*T
      HA(9) = HA(8)*T
      HA(10) = 2.0D0*T*T*X1*X2*X3
      HA(11) = HA(8)
      HA(12) = HA(11)*T
      HA(13) = HA(12)*T
      HA(14) = HA(10)*T
      HA(15) = HA(14)*T
      RETURN

  330 T = 0.1D0*DBLE(KA-1)
      DO 340 I = 1,N* (N+1)/2
          HA(I) = 0.0D0
  340 CONTINUE
      X1 = EXP(-X(2)*T)
      X2 = COS(X(3)*T+X(4))
      X3 = SIN(X(3)*T+X(4))
      X4 = EXP(-X(6)*T)
      HA(2) = -X1*X2*T
      HA(10) = -X1*X2*X(1)
      HA(9) = HA(10)*T
      HA(3) = -HA(9)*T
      HA(7) = -X1*X3
      HA(4) = HA(7)*T
      HA(8) = X1*X3*X(1)*T
      HA(5) = HA(8)*T
      HA(6) = -HA(3)
      HA(20) = -X4*T
      HA(21) = -HA(20)*X(5)*T
      RETURN

  350 BETA = 0.5D0*PI*Y(KA)
      DO 360 I = 1,3
          J = I + I
          XA(I) = X(J-1)
          XB(I) = X(J)
  360 CONTINUE
      CA(4) = CMPLX(1.0D0,0.0D0)
      CB(4) = 1.0D1*CA(4)
      DO 370 J = 1,3
          I = 4 - J
          XC(I) = BETA*XA(I)
          T = XC(I)
          X1 = COS(T)
          X2 = SIN(T)
          C1 = CMPLX(X1,0.0D0)
          C2 = CMPLX(0.0D0, (X2*XB(I)))
          C3 = CMPLX(0.0D0, (X2/XB(I)))
          CB(I) = C1*CB(I+1) + C2*CA(I+1)
          CA(I) = C3*CB(I+1) + C1*CA(I+1)
  370 CONTINUE
      C1 = -CA(1)
      C2 = CB(1) - C1
      C3 = 1.0D0 + 2.0D0*C1/C2
      FA = CDABS(C3)
      C3 = CONJG(C3)
      C1 = 2.0D0/ (C2*C2)
      T = BETA
      DO 380 I = 1,3
          ST(I) = SIN(XC(I))
          CT(I) = COS(XC(I))
          CC(I+I) = (CB(I)*CA(I)-CB(I+1)*CA(I+1))/XB(I)
          CC(I+I-1) = (CB(I)*CA(I+1)-CB(I+1)*CA(I))*T/ST(I)
  380 CONTINUE
      DO 390 I = 1,6
          GA(I) = DBLE(C3*C1*CC(I))/FA
  390 CONTINUE
      CI = CMPLX(0D0,1D0)
      X2 = X(2)
      X4 = X(4)
      X3 = X(6)
      C2 = - (C1+C1)/C2
      S1 = CMPLX(-ST(1)*X2,CT(1))
      S2 = CMPLX(ST(1)/X2,-CT(1))
      S3 = S1*CT(2) - CI*S2*ST(2)*X4
      S4 = S2*CT(2) - CI*S1*ST(2)/X4
      DD(1) = C2* (T*CI* (CB(1)/X2+CA(1)*X2))
      DD(2) = C2* (-ST(1)*CI* (CB(2)/ (X2*X2)-CA(2)))
      DD(3) = C2* (T* (CB(2)*S1/X4-CA(2)*S2*X4))
      DD(4) = C2* (-ST(2)* (CB(3)*S1/ (X4*X4)+CA(3)*S2))
      DD(5) = C2* (T* (CB(3)*S3/X3-CA(3)*S4*X3))
      DD(6) = C2* (-ST(3)* (1D1*S3/ (X3*X3)+S4))
      L = 0
      DO 410 I = 1,6
          L = L + I - 1
          DO 400 J = 1,I
              HA(L+J) = (DBLE(CC(I)* (C3*DD(J)+C1*CONJG(C1*CC(J))))-
     +                  GA(I)*GA(J))/FA
  400     CONTINUE
  410 CONTINUE
      DO 420 I = 1,3
          J = I* (I+I+1)
          HA(J-1) = HA(J-1) + DBLE(C3*C1*CI*T*
     +              (CA(I)*CA(I)+ (CB(I)/XB(I))**2))/FA
          HA(J) = HA(J) + DBLE(C3*C1/XB(I)*
     +            (CI*ST(I)* (CA(I)*CA(I+1)-CB(I)*CB(I+
     +            1)/ (XB(I)*XB(I)))-CC(I+I)))/FA
  420 CONTINUE
      RETURN

  430 T = Y(41+KA)
      BETA = Y(82+KA)
      X1 = (X(1)+ (1.0D0+X(2))*T)**2 + ((1.0D0-X(2))*BETA)**2
      X2 = (X(3)+ (1.0D0+X(4))*T)**2 + ((1.0D0-X(4))*BETA)**2
      X3 = (X(5)+ (1.0D0+X(6))*T)**2 + ((1.0D0-X(6))*BETA)**2
      X4 = (X(7)+ (1.0D0+X(8))*T)**2 + ((1.0D0-X(8))*BETA)**2
      IF (X1.EQ.0.0D0) X1 = 1.0D-30
      IF (X2.EQ.0.0D0) X2 = 1.0D-30
      IF (X3.EQ.0.0D0) X3 = 1.0D-30
      IF (X4.EQ.0.0D0) X4 = 1.0D-30
      FA = SQRT(X1/X2)*SQRT(X3/X4)
      HA(37) = FA/X1* (X(1)+T* (1.0D0+X(2)))
      HA(38) = FA/X1* (X(2)+2.0D0*T*T-1.0D0+X(1)*T)
      HA(39) = -FA/X2* (X(3)+T* (1.0D0+X(4)))
      HA(40) = -FA/X2* (X(4)+2.0D0*T*T-1.0D0+X(3)*T)
      HA(41) = FA/X3* (X(5)+T* (1.0D0+X(6)))
      HA(42) = FA/X3* (X(6)+2.0D0*T*T-1.0D0+X(5)*T)
      HA(43) = -FA/X4* (X(7)+T* (1.0D0+X(8)))
      HA(44) = -FA/X4* (X(8)+2.0D0*T*T-1.0D0+X(7)*T)
      HA(45) = 0.0D0
      FA = X(9)*FA
      DO 440 I = 1,8
          GA(I) = HA(36+I)*X(9)
  440 CONTINUE
      DO 460 J = 1,8
          DO 450 I = 1,J
              HA((J-1)*J/2+I) = GA(I)*GA(J)/FA
  450     CONTINUE
  460 CONTINUE
      HA(1) = FA/X1 - HA(1)
      HA(2) = FA*T/X1 - HA(2)
      HA(3) = FA/X1 - HA(3)
      HA(6) = 3.0D0*HA(6) - FA/X2
      HA(9) = 3.0D0*HA(9) - FA*T/X2
      HA(10) = 3.0D0*HA(10) - FA/X2
      HA(15) = FA/X3 - HA(15)
      HA(20) = FA*T/X3 - HA(20)
      HA(21) = FA/X3 - HA(21)
      HA(28) = 3.0D0*HA(28) - FA/X4
      HA(35) = 3.0D0*HA(35) - FA*T/X4
      HA(36) = 3.0D0*HA(36) - FA/X4
      RETURN

  470 DO 480 I = 1,N* (N+1)/2
          HA(I) = 0.0D0
  480 CONTINUE
      HA(1) = 2.0D0
      HA(3) = 1.0D1
      HA(6) = 1.2D1*X(3)**2
      HA(10) = 6.0D0
      HA(15) = 3.0D2*X(5)**4
      HA(21) = 1.4D1
      HA(27) = -4.0D0
      HA(28) = 1.2D1*X(7)**2
      GO TO (530,490,500,510,520) KA

  490 HA(1) = HA(1) + 4.0D1
      HA(3) = HA(3) + 3.6D2*X(2)**2
      HA(10) = HA(10) + 8.0D1
      RETURN

  500 HA(6) = HA(6) + 2.0D2
      RETURN

  510 HA(3) = HA(3) + 2.0D1
      HA(21) = HA(21) + 1.2D2
      RETURN

  520 HA(1) = HA(1) + 8.0D1
      HA(2) = HA(2) - 3.0D1
      HA(3) = HA(3) + 2.0D1
      HA(6) = HA(6) + 4.0D1
  530 RETURN

  540 DO 550 I = 1,N* (N+1)/2
          HA(I) = 0.0D0
  550 CONTINUE
      HA(1) = 2.0D0
      HA(2) = 1.0D0
      HA(3) = 2.0D0
      HA(6) = 2.0D0
      HA(10) = 8.0D0
      HA(15) = 2.0D0
      HA(21) = 4.0D0
      HA(28) = 1.0D1
      HA(36) = 1.4D1
      HA(45) = 4.0D0
      HA(55) = 2.0D0
      GO TO (610,560,570,580,590,610,610,600,610) KA

  560 HA(1) = HA(1) + 6.0D1
      HA(3) = HA(3) + 8.0D1
      HA(6) = HA(6) + 4.0D1
      RETURN

  570 HA(1) = HA(1) + 1.0D2
      HA(6) = HA(6) + 2.0D1
      RETURN

  580 HA(1) = HA(1) + 1.0D1
      HA(3) = HA(3) + 4.0D1
      HA(15) = HA(15) + 6.0D1
      RETURN

  590 HA(1) = HA(1) + 1.0D1
      HA(2) = HA(2) - 2.0D1
      HA(3) = HA(3) + 4.0D1
      RETURN

  600 HA(45) = HA(45) + 2.4D2
  610 RETURN

  620 DO 630 I = 1,N* (N+1)/2
          HA(I) = 0.0D0
  630 CONTINUE
      HA(1) = 2.0D0
      HA(2) = 1.0D0
      HA(3) = 2.0D0
      HA(6) = 2.0D0
      HA(10) = 8.0D0
      HA(15) = 2.0D0
      HA(21) = 4.0D0
      HA(28) = 1.0D1
      HA(36) = 1.4D1
      HA(45) = 4.0D0
      HA(55) = 2.0D0
      HA(66) = 2.0D0
      HA(78) = 2.0D1
      HA(91) = 1.0D1
      HA(105) = 8.0D0
      HA(120) = 5.4D1
      HA(136) = 1.2D1*X(16)**2
      HA(153) = 2.0D0
      HA(171) = 2.6D1
      HA(190) = 2.0D0
      HA(210) = 2.0D0
      GO TO (610,560,570,580,590,610,610,600,610,610,640,650,660,
     +       670,680,690,700,710) KA

  640 HA(1) = HA(1) + 2.0D1
      RETURN

  650 HA(91) = HA(91) + 1.0D2
      RETURN

  660 HA(91) = HA(91) + 6.0D1
      RETURN

  670 HA(1) = HA(1) + 2.8D2
      RETURN

  680 HA(3) = HA(3) + 3.0D2
      RETURN

  690 HA(1) = HA(1) + 1.0D2
      HA(153) = HA(153) + 10.8D2*X(17)**2
      RETURN

  700 HA(1) = HA(1) + 2.0D1
      RETURN

  710 HA(1) = HA(1) + 1.4D2
      HA(3) = HA(3) + 1.0D2
      HA(190) = HA(190) + 2.0D1
      RETURN

  720 X1 = 0.0D0
      DO 730 I = 1,N
          X3 = 1.0D0
          X4 = X(I)
          IF (I.EQ.1) X3 = 1.0D-8
          IF (I.EQ.4) X3 = 4.0D0
          IF (I.EQ.2 .AND. KA.EQ.1) X4 = X(I) + 2.0D0
          IF (I.EQ.2 .AND. KA.EQ.2) X4 = X(I) - 2.0D0
          X1 = X1 + X3*X4**2
  730 CONTINUE
      X2 = EXP(X1)
      L = 0
      DO 750 I = 1,N
          X3 = 2.0D0
          X4 = X(I)
          IF (I.EQ.1) X3 = 2.0D-8
          IF (I.EQ.4) X3 = 8.0D0
          IF (I.EQ.2 .AND. KA.EQ.1) X4 = X(I) + 2.0D0
          IF (I.EQ.2 .AND. KA.EQ.2) X4 = X(I) - 2.0D0
          DO 740 J = 1,I
              L = L + 1
              X5 = 2.0D0
              X6 = X(J)
              IF (J.EQ.1) X5 = 2.0D-8
              IF (J.EQ.4) X5 = 8.0D0
              IF (J.EQ.2 .AND. KA.EQ.1) X6 = X(J) + 2.0D0
              IF (J.EQ.2 .AND. KA.EQ.2) X6 = X(J) - 2.0D0
              HA(L) = X2*X3*X4*X5*X6
  740     CONTINUE
          HA(L) = HA(L) + X2*X3
  750 CONTINUE
      RETURN

  760 DO 770 I = 1,N* (N+1)/2
          HA(I) = 0.0D0
  770 CONTINUE
      L = 0
      DO 780 I = 1,N
          L = L + I
          X1 = 1.0D0*DBLE(I+KA-1)
          X2 = (X(I)-SIN(DBLE(2*I+KA-3)))**2
          HA(L) = 2.0D0*X1*EXP(X2)* (1.0D0+2.0D0*X2)
  780 CONTINUE
      RETURN

  790 IF (KA.LE.2) THEN
          DO 800 I = 1,N* (N+1)/2
              HA(I) = 0.0D0
  800     CONTINUE
          IF (KA.EQ.1) THEN

          ELSE IF (KA.EQ.2) THEN
              HA(1) = -2.0D0
          END IF

      ELSE
          T = DBLE(KA-2)/2.9D1
          L = 0
          DO 820 I = 1,N
              DO 810 J = 1,I
                  L = L + 1
                  HA(L) = -2.0D0*T** (I+J-2)
  810         CONTINUE
  820     CONTINUE
      END IF

      RETURN

  830 T = 1.0D-1*DBLE(KA-1)
      DO 840 I = 1,N* (N+1)/2
          HA(I) = 0.0D0
  840 CONTINUE
      X1 = EXP(-X(5)*T)
      X2 = EXP(-X(6)* (T-X(9))**2)
      X3 = EXP(-X(7)* (T-X(10))**2)
      X4 = EXP(-X(8)* (T-X(11))**2)
      HA(11) = X1*T
      HA(15) = -X1*X(1)*T**2
      HA(17) = X2* (T-X(9))**2
      HA(21) = -HA(17)*X(2)* (T-X(9))**2
      HA(24) = X3* (T-X(10))**2
      HA(28) = -HA(24)*X(3)* (T-X(10))**2
      HA(32) = X4* (T-X(11))**2
      HA(36) = -HA(32)*X(4)* (T-X(11))**2
      HA(38) = -2.0D0*X2*X(6)* (T-X(9))
      HA(42) = -2.0D0*X2*X(2)* (T-X(9)) +
     +         2.0D0*X2*X(2)*X(6)* (T-X(9))**3
      HA(45) = 2.0D0*X2*X(2)*X(6) - 4.0D0*X2*X(2)* (X(6)* (T-X(9)))**2
      HA(48) = -2.0D0*X3*X(7)* (T-X(10))
      HA(52) = -2.0D0*X3*X(3)* (T-X(10)) +
     +         2.0D0*X3*X(3)*X(7)* (T-X(10))**3
      HA(55) = 2.0D0*X3*X(3)*X(7) - 4.0D0*X3*X(3)* (X(7)* (T-X(10)))**2
      HA(59) = -2.0D0*X4*X(8)* (T-X(11))
      HA(63) = -2.0D0*X4*X(4)* (T-X(11)) +
     +         2.0D0*X4*X(4)*X(8)* (T-X(11))**3
      HA(66) = 2.0D0*X4*X(4)*X(8) - 4.0D0*X4*X(4)* (X(8)* (T-X(11)))**2
      RETURN

      END
* SUBROUTINE TIUD19                ALL SYSTEMS                99/12/01
C PORTABILITY : ALL SYSTEMS
C 94/12/01 VL : ORIGINAL VERSION
*
* PURPOSE :
*  INITIATION OF VARIABLES FOR NONSMOOTH OPTIMIZATION.
*  UNCONSTRAINED DENSE VERSION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  RO  X(N)  VECTOR OF VARIABLES.
*  RO  FMIN  LOWER BOUND FOR VALUE OF THE OBJECTIVE FUNCTION.
*  RO  XMAX  MAXIMUM STEPSIZE.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*  IO  IERR  ERROR INDICATOR.
*
      SUBROUTINE TIUD19(N,X,FMIN,XMAX,NEXT,IERR)
C     .. Parameters ..
      DOUBLE PRECISION ETA9
      PARAMETER (ETA9=1.0D60)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION FMIN,XMAX
      INTEGER IERR,N,NEXT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION X(N)
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION Y(2700)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AI,AJ,AK
      INTEGER I,J,K,KK,L
C     ..
C     .. Local Arrays ..
      INTEGER AA(59),CC(95),PP(23)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,COS,DBLE,EXP,SIN
C     ..
C     .. Common blocks ..
      COMMON /EMPR19/Y
C     ..
C     .. Data statements ..
      DATA AA/5*0,2,3*1,3,1,2,1,1,2,1,4,1,2,2,3,2,1,0,1,0,2,1,0,7*1,0,1,
     +     2,1,0,0,2,1,0,1,1,2,0,0,1,5,10,2,4,3,3*6/
      DATA PP/0,2,3,4,5,6,2,3,-1,4*2,1,1,5,4*1,2,3,2/
      DATA CC/-16,4*0,2,2*-1,2*1,2,-2,0,-2,-9,0,-1,-2,2,1,2*0,2,0,-2,-4,
     +     -1,-3,3,1,1,4,0,-4,1,0,-1,-2,4,1,0,2,0,-1,2*0,2*-1,5,1,-40,
     +     -2,0,2*-4,-1,-40,-60,5,1,30,-20,-10,32,-10,-20,39,-6,-31,32,
     +     -10,-6,10,-6,-10,32,-31,-6,39,-20,-10,32,-10,-20,30,4,8,10,6,
     +     2,-15,-27,-36,-18,-12/
C     ..
      FMIN = 0.0D0
      XMAX = 1.0D3
      IERR = 0
      IF (N.LT.2) GO TO 450
      DO 10 I = 1,N
          X(I) = 0D0
   10 CONTINUE
      GO TO (20,30,40,50,60,100,110,130,140,
     +       150,160,170,190,220,230,240,300,320,
     +       350,350,370,400,420,420,430) NEXT

   20 N = 2
      X(1) = -1.2D0
      X(2) = 1.0D0
      RETURN

   30 N = 2
      X(1) = -1.5D0
      X(2) = 2.0D0
      RETURN

   40 N = 2
      X(1) = 1.0D0
      X(2) = -0.1D0
      RETURN

   50 N = 2
      X(1) = 2.0D0
      X(2) = 2.0D0
      RETURN

   60 N = 2
   70 FMIN = -ETA9
   80 DO 90 I = 1,N
          X(I) = 1.0D0
   90 CONTINUE
      RETURN

  100 N = 2
      X(1) = -1.0D0
      X(2) = 5.0D0
      RETURN

  110 N = 2
      X(1) = -0.5D0
      X(2) = -0.5D0
  120 FMIN = -ETA9
      RETURN

  130 N = 2
      X(1) = 0.8D0
      X(2) = 0.6D0
      GO TO 120

  140 N = 2
      X(1) = -1.0D0
      X(2) = -1.0D0
      GO TO 120

  150 N = 2
      X(1) = 3.0D0
      X(2) = 2.0D0
      GO TO 120

  160 IF (N.LT.4) GO TO 450
      N = 4
      GO TO 120

  170 IF (N.LT.5) GO TO 450
      N = 5
      X(5) = 1.0D0
      DO 180 I = 1,59
          Y(I) = AA(I)
  180 CONTINUE
      Y(57) = 1.7D0
      Y(58) = 2.5D0
      Y(60) = 3.5D0
      RETURN

  190 IF (N.LT.5) GO TO 450
      N = 5
      X(5) = 1.0D0
      FMIN = -ETA9
  200 DO 210 I = 1,95
          Y(I) = CC(I)
  210 CONTINUE
      Y(3) = -3.5D0
      Y(45) = -2.8D0
      Y(53) = -.25D0
      RETURN

  220 IF (N.LT.5) GO TO 450
      N = 5
      X(1) = -2.0D0
      X(2) = 1.5D0
      X(3) = 2.0D0
      X(4) = -1.0D0
      X(5) = -1.0D0
      FMIN = -ETA9
      XMAX = 1.0D0
      RETURN

  230 IF (N.LT.6) GO TO 450
      N = 6
      XMAX = 2D0
      X(1) = 2.0D0
      X(2) = 2.0D0
      X(3) = 7.0D0
      X(5) = -2.0D0
      X(6) = 1.0D0
      RETURN

  240 IF (N.LT.10) GO TO 450
      N = 10
      KK = 0
      DO 290 K = 1,5
          AK = DBLE(K)
          DO 260 I = 1,N
              AI = DBLE(I)
              DO 250 J = I,N
                  AJ = DBLE(J)
                  Y(KK+ (I-1)*N+J) = EXP(AI/AJ)*COS(AI*AJ)*SIN(AK)
                  Y(KK+ (J-1)*N+I) = Y(KK+ (I-1)*N+J)
  250         CONTINUE
  260     CONTINUE
          DO 280 I = 1,N
              AI = DBLE(I)
              Y(100+KK+I) = EXP(AI/AK)*SIN(AI*AK)
              L = KK + (I-1)*N + I
              Y(L) = ABS(SIN(AK))*AI/DBLE(N)
              DO 270 J = 1,N
                  IF (J.NE.I) Y(L) = Y(L) + ABS(Y(KK+ (I-1)*N+J))
  270         CONTINUE
  280     CONTINUE
          KK = KK + 110
  290 CONTINUE
      GO TO 70

  300 IF (N.LT.10) GO TO 450
      N = 10
      XMAX = 1.0D1
      DO 310 I = 1,N
          X(I) = -0.1D0
  310 CONTINUE
      RETURN

  320 IF (N.LT.12) GO TO 450
      N = 12
      DO 330 I = 1,23
          Y(I) = PP(I)
  330 CONTINUE
      Y(10) = -0.5D0
      X(1) = 2.0D0/3.0D0
      X(7) = 5.0D0/3.0D0
      DO 340 I = 2,5
          X(I) = (X(I-1)+Y(I)+Y(I+1))/3.0D0
          X(I+6) = (X(I+5)+Y(I+6)+Y(I+7))/3.0D0
  340 CONTINUE
      X(6) = (X(5)+11.5D0)/3.0D0
      X(12) = (X(11)+1.0D0)/3.0D0
      RETURN

  350 IF (N.LT.20) GO TO 450
      N = 20
      IF (NEXT.EQ.19) XMAX = 5.0D0
      IF (NEXT.EQ.20) XMAX = 1.0D1
      DO 360 I = 1,10
          X(I) = DBLE(I)
          X(I+10) = DBLE(-I-10)
  360 CONTINUE
      RETURN

  370 IF (N.LT.48) GO TO 450
      N = 48
      OPEN (5,FILE='test19.dat')
      READ (5,FMT=*) ((Y(N* (I-2)+J),J=I,N),I=2,N), (Y(N*N+I),I=1,N),
     +  (Y(N*N+N+I),I=1,N)
      DO 390 I = 1,N
          Y(N* (I-1)+I) = 1.0D5
          DO 380 J = 1,I - 1
              Y(N* (I-1)+J) = Y(N* (J-1)+I)
  380     CONTINUE
  390 CONTINUE
      CLOSE (5)
      GO TO 120

  400 IF (N.LT.50) GO TO 450
      N = 50
      XMAX = 1.0D1
      DO 410 I = 1,N
          X(I) = DBLE(I) - 25.5D0
  410 CONTINUE
      RETURN

  420 IF (N.LT.30) GO TO 450
      IF (N.GE.50) N = 50
      IF (N.LT.50) N = 30
      IF (NEXT.EQ.23) XMAX = 5.0D0
      GO TO 80

  430 IF (N.LT.15) GO TO 450
      N = 15
      DO 440 I = 1,N
          X(I) = 1.0D-4
  440 CONTINUE
      X(12) = 6.0D1
      GO TO 200

  450 IERR = 1
      RETURN

      END
* SUBROUTINE TFFU19                ALL SYSTEMS                99/12/01
C PORTABILITY : ALL SYSTEMS
C 94/12/01 RA : ORIGINAL VERSION
*
* PURPOSE :
*  VALUE OF THE NONSMOOTH OBJECTIVE FUNCTION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  F  VALUE OF THE OBJECTIVE FUNCTION.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*
      SUBROUTINE TFFU19(N,X,F,NEXT)
C     .. Parameters ..
      DOUBLE PRECISION ETA9
      PARAMETER (ETA9=1.0D60)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION F
      INTEGER N,NEXT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION X(N)
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION Y(2700)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AI,AJ,F1,F2,F3,F4,T,Z
      INTEGER I,J,K,KK,LL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,COS,DBLE,EXP,MAX,SIN,SQRT
C     ..
C     .. Common blocks ..
      COMMON /EMPR19/Y
C     ..
      F = 0.0D0
      GO TO (10,20,30,40,50,60,70,80,90,
     +       100,110,120,160,210,220,240,300,370,
     +       400,420,440,480,500,540,570) NEXT

   10 F = 100.0D0* (X(2)-X(1)**2)**2 + (1.0D0-X(1))**2
      RETURN

   20 T = X(1)**2 + (X(2)-1.0D0)**2 - 1.0D0
      F = X(2) + ABS(T)
      RETURN

   30 F1 = X(1)**2 + X(2)**4
      F2 = (2.0D0-X(1))**2 + (2.0D0-X(2))**2
      F3 = 2.0D0*EXP(-X(1)+X(2))
      F = MAX(F1,F2,F3)
      RETURN

   40 F1 = X(1)**4 + X(2)**2
      F2 = (2.0D0-X(1))**2 + (2.0D0-X(2))**2
      F3 = 2.0D0*EXP(-X(1)+X(2))
      F = MAX(F1,F2,F3)
      RETURN

   50 F1 = 5.0D0*X(1) + X(2)
      F2 = -5.0D0*X(1) + X(2)
      F3 = X(1)**2 + X(2)**2 + 4.0D0*X(2)
      F = MAX(F1,F2,F3)
      RETURN

   60 F1 = X(1)**2 + X(2)**2
      F2 = F1 + 10.0D0* (-4.0D0*X(1)-X(2)+4.0D0)
      F3 = F1 + 10.0D0* (-X(1)-2.0D0*X(2)+6.0D0)
      F = MAX(F1,F2,F3)
      RETURN

   70 F1 = -X(1) - X(2)
      F2 = F1 + (X(1)**2+X(2)**2-1.0D0)
      F = MAX(F1,F2)
      RETURN

   80 F1 = X(1)**2 + X(2)**2 - 1.0D0
      F2 = 0.0D0
      F = MAX(F1,F2)
      F = 20.0D0*F - X(1)
      RETURN

   90 F1 = X(1)**2 + X(2)**2 - 1.0D0
      IF (F1.LT.0.0D0) F1 = -F1
      F = -X(1) + 2.0D0* (X(1)**2+X(2)**2-1.0D0) + 1.75D0*F1
      RETURN

  100 IF (X(1).GT.ABS(X(2))) THEN
          F = 5.0D0*SQRT(9.0D0*X(1)**2+16.0D0*X(2)**2)

      ELSE IF ((X(1).GT.0.0D0) .AND. (X(1).LE.ABS(X(2)))) THEN
          F = 9.0D0*X(1) + 16.0D0*ABS(X(2))

      ELSE IF (X(1).LE.0.0D0) THEN
          F = 9.0D0*X(1) + 16.0D0*ABS(X(2)) - X(1)**9
      END IF

      RETURN

  110 F = X(1)**2 + X(2)**2 + X(3)**2
      F1 = F + X(3)**2 + X(4)**2 - 5.0D0* (X(1)+X(2)) - 21.0D0*X(3) +
     +     7.0D0*X(4)
      F2 = F + X(4)**2 + X(1) - X(2) + X(3) - X(4) - 8.0D0
      F3 = F + X(2)**2 + 2.0D0*X(4)**2 - X(1) - X(4) - 10.0D0
      F4 = F + 2.0D0*X(1) - X(2) - X(4) - 5.0D0
      F = F1 + 10.0D0*MAX(0.0D0,F2,F3,F4)
      RETURN

  120 F = 0.0D0
      DO 130 J = 1,5
          F = F + (X(J)-Y(J))**2
  130 CONTINUE
      F = F*Y(50+1)
      DO 150 I = 2,10
          F1 = 0.0D0
          DO 140 J = 1,5
              F1 = F1 + (X(J)-Y(5* (I-1)+J))**2
  140     CONTINUE
          F1 = F1*Y(50+I)
          F = MAX(F,F1)
  150 CONTINUE
      RETURN

  160 F = 0.0D0
      DO 180 I = 1,10
          T = Y(50+I)
          DO 170 J = 1,5
              T = T - Y(I+J*10-10)*X(J)
  170     CONTINUE
          F = MAX(F,T)
  180 CONTINUE
      F = F*5D1
      DO 200 J = 1,5
          T = 0D0
          DO 190 I = 1,5
              T = T + Y(55+I+J*5)*X(I)
  190     CONTINUE
          F = F + Y(85+J)*X(J)**3 + Y(90+J)*X(J) + T*X(J)
  200 CONTINUE
      RETURN

  210 F1 = X(1)*X(2)*X(3)*X(4)*X(5)
      F2 = ABS(X(1)**2+X(2)**2+X(3)**2+X(4)**2+X(5)**2-1.0D1)
      F2 = F2 + ABS(X(2)*X(3)-5.0D0*X(4)*X(5))
      F2 = F2 + ABS(X(1)**3+X(2)**3+1.0D0)
      F = F1 + 1.0D1*F2
      RETURN

  220 F = 0.0D0
      DO 230 I = 1,51
          T = DBLE(I-1)/10.0D0
          Z = 0.5D0*EXP(-T) - EXP(-T*2.0D0) + 0.5D0*EXP(-T*3.0D0) +
     +        1.5D0*EXP(-T*1.5D0)*SIN(7.0D0*T) +
     +        EXP(-T*2.5D0)*SIN(5.0D0*T)
          F1 = EXP(-X(2)*T)
          F2 = EXP(-X(6)*T)
          F3 = COS(X(3)*T+X(4))
          F = F + ABS(X(1)*F1*F3+X(5)*F2-Z)
  230 CONTINUE
      RETURN

  240 F = 0.0D0
      K = 1
      DO 260 I = 1,N
          F = F - Y(100+I)*X(I)
          DO 250 J = 1,N
              F = F + Y((I-1)*N+J)*X(I)*X(J)
  250     CONTINUE
  260 CONTINUE
      KK = 110
      DO 290 K = 2,5
          F1 = 0.0D0
          DO 280 I = 1,N
              F1 = F1 - Y(KK+100+I)*X(I)
              DO 270 J = 1,N
                  F1 = F1 + Y(KK+ (I-1)*N+J)*X(I)*X(J)
  270         CONTINUE
  280     CONTINUE
          F = MAX(F,F1)
          KK = KK + 110
  290 CONTINUE
      RETURN

  300 F1 = 0.0D0
      DO 310 I = 1,10
          F1 = F1 + X(I)**2
  310 CONTINUE
      F4 = F1 - 0.25D0
      F1 = 1.0D-3*F4*F4
      DO 320 I = 1,10
          F1 = F1 + (X(I)-1.0D0)**2
  320 CONTINUE
      F2 = 0.0D0
      DO 350 I = 2,30
          AI = DBLE(I-1)/29D0
          F = 0.0D0
          DO 330 J = 1,10
              F = F + X(J)*AI** (J-1)
  330     CONTINUE
          F = -F*F - 1.0D0
          DO 340 J = 2,10
              AJ = DBLE(J-1)
              F = F + X(J)*AJ*AI** (J-2)
  340     CONTINUE
          F2 = F2 + F*F
  350 CONTINUE
      F2 = F2 + X(1)**2 + (X(2)-X(1)**2-1.0D0)**2
      F3 = 0.0D0
      DO 360 I = 2,10
          F3 = F3 + 100.0D0* (X(I)-X(I-1)**2)**2 + (1.0D0-X(I))**2
  360 CONTINUE
      F = MAX(F1,F2,F3)
      RETURN

  370 F = SQRT(X(1)**2+X(7)**2) + SQRT((5.5D0-X(6))**2+
     +    (1.0D0+X(12))**2)
      DO 380 J = 1,6
          F = F + Y(12+J)*SQRT((Y(J)-X(J))**2+ (Y(6+J)-X(6+J))**2)
  380 CONTINUE
      DO 390 J = 1,5
          F = F + Y(18+J)*SQRT((X(J)-X(J+1))**2+ (X(J+6)-X(J+7))**2)
  390 CONTINUE
      RETURN

  400 F = X(1)**2
      DO 410 I = 2,20
          F1 = X(I)**2
          F = MAX(F,F1)
  410 CONTINUE
      RETURN

  420 F = ABS(X(1))
      DO 430 I = 2,20
          F1 = ABS(X(I))
          F = MAX(F,F1)
  430 CONTINUE
      RETURN

  440 N = 48
      F = 0.0D0
      KK = N*N
      LL = N*N + N
      DO 450 I = 1,N
          F = F + Y(LL+I)*X(I)
  450 CONTINUE
      DO 470 J = 1,N
          Z = ETA9
          DO 460 I = 1,N
              T = Y((I-1)*N+J) - X(I)
              IF (T.GE.Z) GO TO 460
              Z = T
              K = I
  460     CONTINUE
          F = F + Y(KK+J)*Z
  470 CONTINUE
      F = -F
      RETURN

  480 F1 = 0.0D0
      F = -ETA9
      DO 490 I = 1,50
          F2 = X(I)
          F1 = F1 + F2
          F = MAX(F,F2)
  490 CONTINUE
      F = 5.0D1*F - F1
      RETURN

  500 F = 0.0D0
      DO 510 J = 1,N
          F = F + X(J)/DBLE(1+J-1)
  510 CONTINUE
      F = ABS(F)
      DO 530 I = 2,N
          F1 = 0.0D0
          DO 520 J = 1,N
              F1 = F1 + X(J)/DBLE(I+J-1)
  520     CONTINUE
          F1 = ABS(F1)
          F = MAX(F1,F)
  530 CONTINUE
      RETURN

  540 F = 0.0D0
      DO 560 J = 1,N
          F1 = 0.0D0
          DO 550 I = 1,N
              F1 = F1 + X(I)/DBLE(I+J-1)
  550     CONTINUE
          F1 = ABS(F1)
          F = F + F1
  560 CONTINUE
      RETURN

  570 F1 = 0.0D0
      DO 580 J = 1,5
          F1 = F1 + Y(85+J)*X(J)**3
  580 CONTINUE
      F = ABS(F1+F1)
      DO 600 J = 1,5
          T = 0.0D0
          DO 590 I = 1,5
              T = T + Y(55+I+J*5)*X(I)
  590     CONTINUE
          F = F + T*X(J)
  600 CONTINUE
      DO 610 J = 6,15
          F = F - Y(45+J)*X(J)
  610 CONTINUE
      DO 630 J = 1,5
          T = -3.0D0*Y(85+J)*X(J)*X(J) - Y(90+J)
          DO 620 I = 1,15
              IF (I.LE.5) T = T - 2.0D0*Y(55+I+J*5)*X(I)
              IF (I.GT.5) T = T + Y(I+J*10-15)*X(I)
  620     CONTINUE
          IF (T.GT.0.0D0) F = F + 1.0D2*T
  630 CONTINUE
      DO 640 I = 1,15
          IF (X(I).LT.0.0D0) F = F - 1.0D2*X(I)
  640 CONTINUE
      RETURN

      END
* SUBROUTINE TFGU19                ALL SYSTEMS                99/12/01
C PORTABILITY : ALL SYSTEMS
C 94/12/01 VL : ORIGINAL VERSION
*
* PURPOSE :
*  GRADIENT OF THE NONSMOOTH OBJECTIVE FUNCTION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  G(N)  GRADIENT OF THE OBJECTIVE FUNCTION.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*
      SUBROUTINE TFGU19(N,X,G,NEXT)
C     .. Parameters ..
      DOUBLE PRECISION ETA9
      PARAMETER (ETA9=1.0D60)
C     ..
C     .. Scalar Arguments ..
      INTEGER N,NEXT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION G(N),X(N)
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION Y(2700)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AI,AJ,F,F1,F2,F3,F4,T
      INTEGER I,J,K,KK,L
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,COS,DBLE,EXP,MAX,SIGN,SIN,SQRT
C     ..
C     .. Common blocks ..
      COMMON /EMPR19/Y
C     ..
      DO 10 I = 1,N
          G(I) = 0.0D0
   10 CONTINUE
      GO TO (20,30,40,40,50,60,70,80,90,
     +       100,110,170,210,260,280,300,360,490,
     +       520,540,560,600,620,660,700) NEXT

   20 G(2) = 200.0D0* (X(2)-X(1)**2)
      G(1) = 2.0D0*X(1)* (1.0D0-G(2)) - 2.0D0
      RETURN

   30 T = X(1)**2 + (X(2)-1.0D0)**2 - 1.0D0
      F = SIGN(2.0D0,T)
      G(1) = F*X(1)
      G(2) = F* (X(2)-1.0D0) + 1.0D0
      RETURN

   40 I = NEXT - 2
      J = 5 - NEXT
      F1 = X(I)**2 + X(J)**4
      F2 = (2.0D0-X(1))**2 + (2.0D0-X(2))**2
      F3 = 2.0D0*EXP(-X(1)+X(2))
      IF ((F1.GE.F2) .AND. (F1.GE.F3)) THEN
          G(I) = 2.0D0*X(I)
          G(J) = 4.0D0*X(J)**3

      ELSE IF ((F2.GE.F1) .AND. (F2.GE.F3)) THEN
          G(1) = -2.0D0* (2.0D0-X(1))
          G(2) = -2.0D0* (2.0D0-X(2))

      ELSE
          G(2) = 2.0D0*EXP(-X(1)+X(2))
          G(1) = -G(2)
      END IF

      RETURN

   50 F1 = 5.0D0*X(1) + X(2)
      F2 = -5.0D0*X(1) + X(2)
      F3 = X(1)**2 + X(2)**2 + 4.0D0*X(2)
      G(2) = 1.0D0
      IF ((F1.GE.F2) .AND. (F1.GE.F3)) THEN
          G(1) = 5.0D0

      ELSE IF ((F2.GE.F1) .AND. (F2.GE.F3)) THEN
          G(1) = -5.0D0

      ELSE
          G(1) = 2.0D0*X(1)
          G(2) = 2.0D0*X(2) + 4.0D0
      END IF

      RETURN

   60 F1 = X(1)**2 + X(2)**2
      F2 = F1 + 10.0D0* (-4.0D0*X(1)-X(2)+4.0D0)
      F3 = F1 + 10.0D0* (-X(1)-2.0D0*X(2)+6.0D0)
      G(1) = 2.0D0*X(1)
      G(2) = 2.0D0*X(2)
      IF ((F1.GE.F2) .AND. (F1.GE.F3)) THEN

      ELSE IF ((F2.GE.F1) .AND. (F2.GE.F3)) THEN
          G(1) = G(1) - 40.0D0
          G(2) = G(2) - 10.0D0

      ELSE
          G(1) = G(1) - 10.0D0
          G(2) = G(2) - 20.0D0
      END IF

      RETURN

   70 F1 = -X(1) - X(2)
      F2 = F1 + (X(1)**2+X(2)**2-1.0D0)
      IF (F1.GE.F2) THEN
          G(1) = -1.0D0
          G(2) = -1.0D0

      ELSE
          G(1) = -1.0D0 + 2.0D0*X(1)
          G(2) = -1.0D0 + 2.0D0*X(2)
      END IF

      RETURN

   80 F1 = X(1)**2 + X(2)**2 - 1.0D0
      G(1) = -1.0D0
      IF (F1.GE.0.0D0) THEN
          G(1) = 40.0D0*X(1) - 1.0D0
          G(2) = 40.0D0*X(2)
      END IF

      RETURN

   90 F1 = SIGN(3.5D0,X(1)**2+X(2)**2-1.0D0) + 4.0D0
      G(1) = F1*X(1) - 1.0D0
      G(2) = F1*X(2)
      RETURN

  100 IF (X(1).GT.ABS(X(2))) THEN
          G(1) = 45.0D0*X(1)/SQRT(9.0D0*X(1)**2+16.0D0*X(2)**2)
          G(2) = 80.0D0*X(2)/SQRT(9.0D0*X(1)**2+16.0D0*X(2)**2)

      ELSE
          G(1) = 9.0D0
          IF (X(1).LT.0.0D0) G(1) = 9.0D0 - 9.0D0*X(1)**8
          G(2) = SIGN(16.0D0,X(2))
      END IF

      RETURN

  110 F = X(1)**2 + X(2)**2 + X(3)**2
      F2 = F + X(4)**2 + X(1) - X(2) + X(3) - X(4) - 8.0D0
      F3 = F + X(2)**2 + 2.0D0*X(4)**2 - X(1) - X(4) - 10.0D0
      F4 = F + 2.0D0*X(1) - X(2) - X(4) - 5.0D0
      L = 1
      IF (F2.GT.0.0D0) L = 2
      IF (F3.GT.MAX(F2,0.0D0)) L = 3
      IF (F4.GT.MAX(F2,F3,0.0D0)) L = 4
      GO TO (120,130,140,150) L

  120 G(1) = 2.0D0*X(1) - 5.0D0
      G(2) = 2.0D0*X(2) - 5.0D0
      G(3) = 4.0D0*X(3) - 21.0D0
      G(4) = 2.0D0*X(4) + 7.0D0
      RETURN

  130 G(1) = 22.0D0*X(1) + 5.0D0
      G(2) = 22.0D0*X(2) - 15.0D0
      G(3) = 24.0D0*X(3) - 11.0D0
      G(4) = 22.0D0*X(4) - 3.0D0
      RETURN

  140 G(1) = 22.0D0*X(1) - 15.0D0
      G(2) = 42.0D0*X(2) - 5.0D0
      G(4) = 42.0D0*X(4) - 3.0D0
      GO TO 160

  150 G(1) = 22.0D0*X(1) + 15.0D0
      G(2) = 22.0D0*X(2) - 15.0D0
      G(4) = 2.0D0*X(4) - 3.0D0
  160 G(3) = 24.0D0*X(3) - 21.0D0
      RETURN

  170 F = -ETA9
      DO 190 I = 1,10
          F1 = 0.0D0
          DO 180 J = 1,5
              F1 = F1 + (X(J)-Y(5* (I-1)+J))**2
  180     CONTINUE
          F1 = F1*Y(50+I)
          IF (F.LT.F1) K = I
          F = MAX(F,F1)
  190 CONTINUE
      DO 200 J = 1,5
          G(J) = 2.0D0*Y(50+K)* (X(J)-Y(5* (K-1)+J))
  200 CONTINUE
      RETURN

  210 F1 = 0.0D0
      DO 230 I = 1,10
          T = Y(50+I)
          DO 220 J = 1,5
              T = T - Y(I+J*10-10)*X(J)
  220     CONTINUE
          IF (T.GT.F1) K = I
          F1 = MAX(F1,T)
  230 CONTINUE
      DO 250 J = 1,5
          T = 0.0D0
          DO 240 I = 1,5
              T = T + Y(55+I+J*5)*X(I)
  240     CONTINUE
          G(J) = 3.0D0*Y(85+J)*X(J)*X(J) + T + T + Y(90+J)
          IF (F1.GT.0.0D0) G(J) = G(J) - 5.0D1*Y(K+J*10-10)
  250 CONTINUE
      RETURN

  260 G(1) = X(2)*X(3)*X(4)*X(5)
      G(2) = X(1)*X(3)*X(4)*X(5)
      G(3) = X(1)*X(2)*X(4)*X(5)
      G(4) = X(1)*X(2)*X(3)*X(5)
      G(5) = X(1)*X(2)*X(3)*X(4)
      F1 = X(1)**2 + X(2)**2 + X(3)**2 + X(4)**2 + X(5)**2 - 1.0D1
      F4 = 1.0D0
      IF (F1.LT.0.0D0) F4 = -F4
      DO 270 I = 1,5
          G(I) = G(I) + 2.0D1*F4*X(I)
  270 CONTINUE
      F2 = X(2)*X(3) - 5.0D0*X(4)*X(5)
      F4 = 1.0D0
      IF (F2.LT.0.0D0) F4 = -F4
      G(2) = G(2) + 1.0D1*F4*X(3)
      G(3) = G(3) + 1.0D1*F4*X(2)
      G(4) = G(4) - 5.0D1*F4*X(5)
      G(5) = G(5) - 5.0D1*F4*X(4)
      F3 = X(1)**3 + X(2)**3 + 1.0D0
      F4 = 1.0D0
      IF (F3.LT.0.0D0) F4 = -F4
      G(1) = G(1) + 3.0D1*F4*X(1)**2
      G(2) = G(2) + 3.0D1*F4*X(2)**2
      RETURN

  280 DO 290 I = 1,51
          T = DBLE(I-1)/10.0D0
          F = 0.5D0*EXP(-T) - EXP(-T*2.0D0) + 0.5D0*EXP(-T*3.0D0) +
     +        1.5D0*EXP(-T*1.5D0)*SIN(7.0D0*T) +
     +        EXP(-T*2.5D0)*SIN(5.0D0*T)
          F1 = EXP(-X(2)*T)
          F2 = EXP(-X(6)*T)
          F3 = COS(X(3)*T+X(4))
          F4 = SIN(X(3)*T+X(4))
          AI = SIGN(1.0D0,X(1)*F1*F3+X(5)*F2-F)
          G(1) = G(1) + AI*F1*F3
          G(2) = G(2) - AI*F1*F3*X(1)*T
          G(3) = G(3) - AI*F1*F4*X(1)*T
          G(4) = G(4) - AI*F1*F4*X(1)
          G(5) = G(5) + AI*F2
          G(6) = G(6) - AI*F2*X(5)*T
  290 CONTINUE
      RETURN

  300 F = -ETA9
      L = 1
      KK = 0
      DO 330 K = 1,5
          F1 = 0.0D0
          DO 320 I = 1,N
              F1 = F1 - Y(KK+100+I)*X(I)
              DO 310 J = 1,N
                  F1 = F1 + Y(KK+ (I-1)*N+J)*X(I)*X(J)
  310         CONTINUE
  320     CONTINUE
          IF (F.LT.F1) L = K
          F = MAX(F,F1)
          KK = KK + 110
  330 CONTINUE
      DO 350 I = 1,N
          G(I) = -Y((L-1)*110+100+I)
          DO 340 J = 1,N
              G(I) = G(I) + 2.0D0*Y((L-1)*110+ (I-1)*N+J)*X(J)
  340     CONTINUE
  350 CONTINUE
      RETURN

  360 F1 = 0.0D0
      DO 370 I = 1,10
          F1 = F1 + X(I)**2
  370 CONTINUE
      F4 = F1 - 0.25D0
      F1 = 1.0D-3*F4*F4
      DO 380 I = 1,10
          F1 = F1 + (X(I)-1.0D0)**2
  380 CONTINUE
      F2 = 0.0D0
      DO 410 I = 2,30
          AI = DBLE(I-1)/29D0
          F = 0.0D0
          DO 390 J = 1,10
              F = F + X(J)*AI** (J-1)
  390     CONTINUE
          F = -F*F - 1.0D0
          DO 400 J = 2,10
              AJ = DBLE(J-1)
              F = F + X(J)*AJ*AI** (J-2)
  400     CONTINUE
          F2 = F2 + F*F
  410 CONTINUE
      F2 = F2 + X(1)**2 + (X(2)-X(1)**2-1.0D0)**2
      F3 = 0.0D0
      DO 420 I = 2,10
          F3 = F3 + 100.0D0* (X(I)-X(I-1)**2)**2 + (1.0D0-X(I))**2
  420 CONTINUE
      IF ((F1.GE.F2) .AND. (F1.GE.F3)) THEN
          DO 430 I = 1,10
              G(I) = 2.0D0*X(I) - 2.0D0 + 4.0D-3*X(I)*F4
  430     CONTINUE

      ELSE IF ((F2.GE.F1) .AND. (F2.GE.F3)) THEN
          DO 470 J = 1,10
              DO 460 I = 2,30
                  AI = DBLE(I-1)/29D0
                  F = 0.0D0
                  DO 440 K = 1,10
                      F = F - X(K)*AI** (K-1)
  440             CONTINUE
                  T = 2.0D0*F*AI** (J-1)
                  IF (J.GE.2) T = T + (J-1)*AI** (J-2)
                  F = -F*F - 1.0D0
                  DO 450 K = 2,10
                      F = F + X(K)* (K-1)*AI** (K-2)
  450             CONTINUE
                  G(J) = G(J) + 2.0D0*F*T
  460         CONTINUE
  470     CONTINUE
          G(1) = G(1) + 2.0D0*X(1) - 4.0D0*X(1)* (X(2)-X(1)**2-1.0D0)
          G(2) = G(2) + 2.0D0* (X(2)-X(1)**2-1.0D0)

      ELSE
          DO 480 I = 1,10
              G(I) = 0.0D0
              IF (I.GE.2) G(I) = G(I) + 2.0D2* (X(I)-X(I-1)**2) -
     +                           2.0D0* (1.0D0-X(I))
              IF (I.LE.9) G(I) = G(I) - 4.0D2*X(I)* (X(I+1)-X(I)**2)
  480     CONTINUE
      END IF

      RETURN

  490 G(1) = X(1)/SQRT(X(1)**2+X(7)**2)
      G(7) = X(7)/SQRT(X(1)**2+X(7)**2)
      T = SQRT((5.5D0-X(6))**2+ (1.0D0+X(12))**2)
      G(6) = - (5.5D0-X(6))/T
      G(12) = (1.0D0+X(12))/T
      DO 500 J = 1,6
          T = SQRT((Y(J)-X(J))**2+ (Y(6+J)-X(6+J))**2)
          G(J) = G(J) - Y(12+J)* (Y(J)-X(J))/T
          G(6+J) = G(6+J) - Y(12+J)* (Y(J+6)-X(J+6))/T
  500 CONTINUE
      DO 510 J = 1,5
          T = SQRT((X(J)-X(J+1))**2+ (X(J+6)-X(J+7))**2)
          G(J) = G(J) + Y(18+J)* (X(J)-X(J+1))/T
          G(J+1) = G(J+1) - Y(18+J)* (X(J)-X(J+1))/T
          G(J+6) = G(J+6) + Y(18+J)* (X(J+6)-X(J+7))/T
          G(J+7) = G(J+7) - Y(18+J)* (X(J+6)-X(J+7))/T
  510 CONTINUE
      RETURN

  520 F = X(1)**2
      K = 1
      DO 530 I = 2,20
          F1 = X(I)**2
          IF (F.LT.F1) K = I
          F = MAX(F,F1)
  530 CONTINUE
      G(K) = 2.0D0*X(K)
      RETURN

  540 F = ABS(X(1))
      K = 1
      DO 550 I = 2,20
          F1 = ABS(X(I))
          IF (F.LT.F1) K = I
          F = MAX(F,F1)
  550 CONTINUE
      G(K) = SIGN(1.0D0,X(K))
      RETURN

  560 KK = N*N
      DO 570 I = 1,N
          G(I) = -Y(KK+N+I)
  570 CONTINUE
      DO 590 J = 1,N
          F = ETA9
          DO 580 I = 1,N
              T = Y((I-1)*N+J) - X(I)
              IF (T.GE.F) GO TO 580
              F = T
              K = I
  580     CONTINUE
          G(K) = G(K) + Y(KK+J)
  590 CONTINUE
      RETURN

  600 F = -ETA9
      DO 610 I = 1,50
          F2 = X(I)
          IF (F.LT.F2) K = I
          F = MAX(F,F2)
          G(I) = -1.0D0
  610 CONTINUE
      G(K) = G(K) + 5D1
      RETURN

  620 F1 = -ETA9
      DO 640 I = 1,N
          F = 0.0D0
          DO 630 J = 1,N
              F = F + X(J)/DBLE(I+J-1)
  630     CONTINUE
          IF (F1.GE.ABS(F)) GO TO 640
          K = I
          AI = SIGN(1D0,F)
          F1 = ABS(F)
  640 CONTINUE
      DO 650 J = 1,N
          G(J) = AI/DBLE(K+J-1)
  650 CONTINUE
      RETURN

  660 DO 690 J = 1,N
          F1 = 0.0D0
          DO 670 I = 1,N
              F1 = F1 + X(I)/DBLE(I+J-1)
  670     CONTINUE
          AJ = SIGN(1.0D0,F1)
          DO 680 I = 1,N
              G(I) = G(I) + AJ/DBLE(I+J-1)
  680     CONTINUE
  690 CONTINUE
      RETURN

  700 F1 = 0.0D0
      DO 710 J = 1,5
          F1 = F1 + Y(85+J)*X(J)**3
  710 CONTINUE
      DO 730 J = 1,5
          T = 0.0D0
          DO 720 I = 1,5
              T = T + Y(55+I+J*5)*X(I)
  720     CONTINUE
          G(J) = SIGN(6D0,F1)*Y(85+J)*X(J)*X(J) + T + T
  730 CONTINUE
      DO 740 J = 6,15
          G(J) = -Y(45+J)
  740 CONTINUE
      DO 770 J = 1,5
          T = -3.0D0*Y(85+J)*X(J)*X(J) - Y(90+J)
          DO 750 I = 1,15
              IF (I.LE.5) T = T - 2.0D0*Y(55+I+J*5)*X(I)
              IF (I.GT.5) T = T + Y(I+J*10-15)*X(I)
  750     CONTINUE
          IF (T.LE.0.0D0) GO TO 770
          G(J) = G(J) - 6.0D2*Y(85+J)*X(J)
          DO 760 I = 1,15
              IF (I.LE.5) G(I) = G(I) - 2.0D2*Y(55+I+J*5)
              IF (I.GT.5) G(I) = G(I) + 1.0D2*Y(I+J*10-15)
  760     CONTINUE
  770 CONTINUE
      DO 780 I = 1,15
          IF (X(I).LT.0.0D0) G(I) = G(I) - 1.0D2
  780 CONTINUE
      RETURN

      END
* SUBROUTINE TFHD19                ALL SYSTEMS                99/12/01
C PORTABILITY : ALL SYSTEMS
C 95/12/01 VL : ORIGINAL VERSION
*
* PURPOSE :
*  HESSIAN MATRIX OF THE NONSMOOTH OBJECTIVE FUNCTION.
*  DENSE VERSION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  H(N*(N+1)/2)  HESSIAN MATRIX OF THE OBJECTIVE FUNCTION.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*
      SUBROUTINE TFHD19(N,X,H,NEXT)
C     .. Parameters ..
      DOUBLE PRECISION ETA9
      PARAMETER (ETA9=1.0D60)
C     ..
C     .. Scalar Arguments ..
      INTEGER N,NEXT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION H(N* (N+1)/2),X(N)
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION Y(2700)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AI,AJ,F,F1,F2,F3,F4,T
      INTEGER I,J,K,KK,L
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,COS,DBLE,EXP,MAX,SIGN,SIN
C     ..
C     .. Common blocks ..
      COMMON /EMPR19/Y
C     ..
C     .. Statement Functions ..
      INTEGER IN
C     ..
C     .. Statement Function definitions ..
      IN(I,J) = (J-1)*J/2 + I
C     ..
      DO 10 I = 1,N* (N+1)/2
          H(I) = 0.0D0
   10 CONTINUE
      GO TO (20,30,40,40,50,60,70,80,90,
     +       100,110,180,220,260,280,300,360,510,
     +       540,560,570,580,590,600,610) NEXT

   20 H(1) = 1200D0*X(1)**2 - 400D0*X(2) + 2.0D0
      H(2) = -400D0*X(1)
      H(3) = 200D0
      RETURN

   30 T = X(1)**2 + (X(2)-1.0D0)**2 - 1.0D0
      F = SIGN(2.0D0,T)
      H(1) = F
      H(3) = F
      RETURN

   40 I = NEXT - 2
      J = 5 - NEXT
      F1 = X(I)**2 + X(J)**4
      F2 = (2.0D0-X(1))**2 + (2.0D0-X(2))**2
      F3 = 2.0D0*EXP(-X(1)+X(2))
      IF ((F1.GE.F2) .AND. (F1.GE.F3)) THEN
          H(NEXT*2-5) = 2.0D0
          H(9-NEXT*2) = 12.0D0*X(J)**2

      ELSE IF ((F2.GE.F1) .AND. (F2.GE.F3)) THEN
          H(1) = 2.0D0
          H(3) = 2.0D0

      ELSE
          H(1) = 2.0D0*EXP(X(2)-X(1))
          H(2) = -H(1)
          H(3) = H(1)
      END IF

      RETURN

   50 F1 = 5.0D0*X(1) + X(2)
      F2 = -5.0D0*X(1) + X(2)
      F3 = X(1)**2 + X(2)**2 + 4.0D0*X(2)
      IF ((F1.GE.F2) .AND. (F1.GE.F3)) THEN

      ELSE IF ((F2.GE.F1) .AND. (F2.GE.F3)) THEN

      ELSE
          H(1) = 2.0D0
          H(3) = 2.0D0
      END IF

      RETURN

   60 H(1) = 2.0D0
      H(3) = 2.0D0
      RETURN

   70 F1 = -X(1) - X(2)
      F2 = F1 + (X(1)**2+X(2)**2-1.0D0)
      IF (F1.LT.F2) THEN
          H(1) = 2.0D0
          H(3) = 2.0D0
      END IF

      RETURN

   80 F1 = X(1)**2 + X(2)**2 - 1.0D0
      IF (F1.GE.0.0D0) THEN
          H(1) = 40.0D0
          H(3) = 40.0D0
      END IF

      RETURN

   90 F1 = SIGN(3.5D0,X(1)**2+X(2)**2-1.0D0) + 4.0D0
      H(1) = F1
      H(3) = F1
      RETURN

  100 IF (X(1).GT.ABS(X(2))) THEN
          F1 = 720D0* (9.0D0*X(1)**2+16.0D0*X(2)**2)** (-1.5D0)
          H(1) = F1*X(2)**2
          H(2) = -F1*X(1)*X(2)
          H(3) = F1*X(1)**2

      ELSE
          IF (X(1).LT.0.0D0) H(1) = -72.0D0*X(1)**8
      END IF

      RETURN

  110 F = X(1)**2 + X(2)**2 + X(3)**2
      F2 = F + X(4)**2 + X(1) - X(2) + X(3) - X(4) - 8.0D0
      F3 = F + X(2)**2 + 2.0D0*X(4)**2 - X(1) - X(4) - 10.0D0
      F4 = F + 2.0D0*X(1) - X(2) - X(4) - 5.0D0
      L = 1
      IF (F2.GT.0.0D0) L = 2
      IF (F3.GT.MAX(F2,0.0D0)) L = 3
      IF (F4.GT.MAX(F2,F3,0.0D0)) L = 4
      GO TO (120,130,140,150) L

  120 H(1) = 2.0D0
      H(3) = 2.0D0
      H(6) = 4.0D0
      H(10) = 2.0D0
      RETURN

  130 H(10) = 22.0D0
      GO TO 160

  140 H(3) = 42.0D0
      H(10) = 42.0D0
      GO TO 170

  150 H(10) = 2.0D0
  160 H(3) = 22.0D0
  170 H(1) = 22.0D0
      H(6) = 24.0D0
      RETURN

  180 F = -ETA9
      DO 200 I = 1,10
          F1 = 0.0D0
          DO 190 J = 1,5
              F1 = F1 + (X(J)-Y(5* (I-1)+J))**2
  190     CONTINUE
          F1 = F1*Y(50+I)
          IF (F.LT.F1) K = I
          F = MAX(F,F1)
  200 CONTINUE
      DO 210 J = 1,5
          H(J* (J+1)/2) = 2.0D0*Y(50+K)
  210 CONTINUE
      RETURN

  220 DO 230 J = 1,5
          H(J* (J+1)/2) = 6.0D0*Y(85+J)*X(J)
  230 CONTINUE
      K = 1
      DO 250 I = 1,N
          DO 240 J = 1,I
              H(K) = H(K) + 2.0D0*Y(55+I+J*5)
              K = K + 1
  240     CONTINUE
  250 CONTINUE
      RETURN

  260 H(1) = 0.0D0
      H(2) = X(3)*X(4)*X(5)
      H(3) = 0.0D0
      H(4) = X(2)*X(4)*X(5)
      H(5) = X(1)*X(4)*X(5)
      H(6) = 0.0D0
      H(7) = X(2)*X(3)*X(5)
      H(8) = X(1)*X(3)*X(5)
      H(9) = X(1)*X(2)*X(5)
      H(10) = 0.0D0
      H(11) = X(2)*X(3)*X(4)
      H(12) = X(1)*X(3)*X(4)
      H(13) = X(1)*X(2)*X(4)
      H(14) = X(1)*X(2)*X(3)
      H(15) = 0.0D0
      F1 = X(1)**2 + X(2)**2 + X(3)**2 + X(4)**2 + X(5)**2 - 1.0D1
      F4 = 1.0D0
      IF (F1.LT.0.0D0) F4 = -F4
      L = 0
      DO 270 I = 1,5
          L = L + I
          H(L) = H(L) + 2.0D1*F4
  270 CONTINUE
      F2 = X(2)*X(3) - 5.0D0*X(4)*X(5)
      F4 = 1.0D0
      IF (F2.LT.0.0D0) F4 = -F4
      H(5) = H(5) + 1.0D1*F4
      H(14) = H(14) - 5.0D1*F4
      F3 = X(1)**3 + X(2)**3 + 1.0D0
      F4 = 1.0D0
      IF (F3.LT.0.0D0) F4 = -F4
      H(1) = H(1) + 6.0D1*F4*X(1)
      H(3) = H(3) + 6.0D1*F4*X(2)
      RETURN

  280 DO 290 I = 1,51
          T = DBLE(I-1)/10.0D0
          F = 0.5D0*EXP(-T) - EXP(-T*2.0D0) + 0.5D0*EXP(-T*3.0D0) +
     +        1.5D0*EXP(-T*1.5D0)*SIN(7.0D0*T) +
     +        EXP(-T*2.5D0)*SIN(5.0D0*T)
          F1 = EXP(-X(2)*T)
          F2 = EXP(-X(6)*T)
          F3 = COS(X(3)*T+X(4))
          F4 = SIN(X(3)*T+X(4))
          AI = SIGN(1.0D0,X(1)*F1*F3+X(5)*F2-F)
          H(2) = H(2) - AI*F1*F3*T
          H(3) = H(3) + AI*F1*F3*T*T*X(1)
          H(4) = H(4) - AI*F1*F4*T
          H(5) = H(5) + AI*F1*F4*T*T*X(1)
          H(6) = H(6) - AI*F1*F3*T*T*X(1)
          H(7) = H(7) - AI*F1*F4
          H(8) = H(8) + AI*F1*F4*T*X(1)
          H(9) = H(9) - AI*F1*F3*T*X(1)
          H(10) = H(10) - AI*F1*F3*X(1)
          H(15) = H(15) - AI*F2*T
          H(20) = H(20) - AI*F2*T
          H(21) = H(21) + AI*F2*T*T*X(5)
  290 CONTINUE
      RETURN

  300 F = -ETA9
      L = 1
      KK = 0
      DO 330 K = 1,5
          F1 = 0.0D0
          DO 320 I = 1,N
              F1 = F1 - Y(KK+100+I)*X(I)
              DO 310 J = 1,N
                  F1 = F1 + Y(KK+ (I-1)*N+J)*X(I)*X(J)
  310         CONTINUE
  320     CONTINUE
          IF (F.LT.F1) L = K
          F = MAX(F,F1)
          KK = KK + 110
  330 CONTINUE
      K = 1
      DO 350 I = 1,N
          DO 340 J = 1,I
              H(K) = 2.0D0*Y((L-1)*110+ (I-1)*N+J)
              K = K + 1
  340     CONTINUE
  350 CONTINUE
      RETURN

  360 F1 = 0.0D0
      DO 370 I = 1,10
          F1 = F1 + X(I)**2
  370 CONTINUE
      F4 = F1 - 0.25D0
      F1 = 1.0D-3*F4*F4
      DO 380 I = 1,10
          F1 = F1 + (X(I)-1.0D0)**2
  380 CONTINUE
      F2 = 0.0D0
      DO 410 I = 2,30
          AI = DBLE(I-1)/29D0
          F = 0.0D0
          DO 390 J = 1,10
              F = F + X(J)*AI** (J-1)
  390     CONTINUE
          F = -F*F - 1.0D0
          DO 400 J = 2,10
              AJ = DBLE(J-1)
              F = F + X(J)*AJ*AI** (J-2)
  400     CONTINUE
          F2 = F2 + F*F
  410 CONTINUE
      F2 = F2 + X(1)**2 + (X(2)-X(1)**2-1.0D0)**2
      F3 = 0.0D0
      DO 420 I = 2,10
          F3 = F3 + 100.0D0* (X(I)-X(I-1)**2)**2 + (1.0D0-X(I))**2
  420 CONTINUE
      IF ((F1.GE.F2) .AND. (F1.GE.F3)) THEN
          L = 1
          DO 440 I = 1,N
              DO 430 J = 1,I
                  H(L) = 8.0D-3*X(I)*X(J)
                  IF (J.EQ.I) H(L) = H(L) + 2.0D0 + 4.0D-3*F4
                  L = L + 1
  430         CONTINUE
  440     CONTINUE

      ELSE IF ((F2.GE.F1) .AND. (F2.GE.F3)) THEN
          KK = 1
          DO 490 J = 1,N
              DO 480 L = 1,J
                  DO 470 I = 2,30
                      AI = DBLE(I-1)/29.0D0
                      F = 0.0D0
                      DO 450 K = 1,10
                          F = F - X(K)*AI** (K-1)
  450                 CONTINUE
                      T = 2.0D0*F*AI** (J-1)
                      IF (J.GE.2) T = T + (J-1)*AI** (J-2)
                      F4 = 2.0D0*F*AI** (L-1)
                      IF (L.GE.2) F4 = F4 + (L-1)*AI** (L-2)
                      F = -F*F - 1.0D0
                      DO 460 K = 2,10
                          F = F + X(K)* (K-1)*AI** (K-2)
  460                 CONTINUE
                      H(KK) = H(KK) + 2.0D0* (T*F4-2.0D0*F*AI** (J+L-2))
  470             CONTINUE
                  KK = KK + 1
  480         CONTINUE
  490     CONTINUE
          H(1) = H(1) + 12.0D0*X(1)**2 - 4.0D0*X(2) + 6.0D0
          H(2) = H(2) - 4.0D0*X(1)
          H(3) = H(3) + 2.0D0

      ELSE
          DO 500 I = 1,10
              J = I* (I+1)/2
              IF (I.GE.2) H(J) = 202.0D0
              IF (I.LE.9) H(J) = H(J) + 4.0D2* (3.0D0*X(I)**2-X(I+1))
              IF (I.LE.9) H(J+I) = -4.0D2*X(I)
  500     CONTINUE
      END IF

      RETURN

  510 T = (X(1)**2+X(7)**2)**1.5D0
      H(1) = X(7)**2/T
      H(22) = -X(1)*X(7)/T
      H(28) = X(1)**2/T
      T = ((5.5D0-X(6))**2+ (1.0D0+X(12))**2)**1.5D0
      H(21) = (1.0D0+X(12))**2/T
      H(72) = (5.5D0-X(6))* (1.0D0+X(12))/T
      H(78) = (5.5D0-X(6))**2/T
      DO 520 J = 1,6
          T = Y(12+J)/ ((Y(J)-X(J))**2+ (Y(6+J)-X(6+J))**2)**1.5D0
          H(IN(J,J)) = H(IN(J,J)) + (Y(J+6)-X(J+6))**2*T
          H(IN(J,J+6)) = H(IN(J,J+6)) - (Y(J+6)-X(J+6))* (Y(J)-X(J))*T
          H(IN(J+6,J+6)) = H(IN(J+6,J+6)) + (Y(J)-X(J))**2*T
  520 CONTINUE
      DO 530 J = 1,6
          IF (J.LT.6) THEN
              F1 = X(J) - X(J+1)
              F2 = X(J+6) - X(J+7)
              T = Y(18+J)/ (F1*F1+F2*F2)**1.5D0
              H(IN(J,J)) = H(IN(J,J)) + F2*F2*T
              H(IN(J,J+1)) = H(IN(J,J+1)) - F2*F2*T
              H(IN(J,J+6)) = H(IN(J,J+6)) - F1*F2*T
              H(IN(J,J+7)) = H(IN(J,J+7)) + F1*F2*T
              H(IN(J+6,J+6)) = H(IN(J+6,J+6)) + F1*F1*T
              H(IN(J+6,J+7)) = H(IN(J+6,J+7)) - F1*F1*T
          END IF

          IF (J.GT.1) THEN
              F1 = X(J) - X(J-1)
              F2 = X(J+6) - X(J+5)
              T = Y(17+J)/ (F1*F1+F2*F2)**1.5D0
              H(IN(J,J)) = H(IN(J,J)) + F2*F2*T
              H(IN(J,J+5)) = H(IN(J,J+5)) + F1*F2*T
              H(IN(J,J+6)) = H(IN(J,J+6)) - F1*F2*T
              H(IN(J+6,J+6)) = H(IN(J+6,J+6)) + F1*F1*T
          END IF

  530 CONTINUE
      RETURN

  540 F = X(1)**2
      K = 1
      DO 550 I = 2,20
          F1 = X(I)**2
          IF (F.LT.F1) K = I
          F = MAX(F,F1)
  550 CONTINUE
      H(K* (K+1)/2) = 2.0D0
      RETURN

  560 RETURN

  570 RETURN

  580 RETURN

  590 RETURN

  600 RETURN

  610 F1 = 0D0
      DO 620 J = 1,5
          F1 = F1 + Y(85+J)*X(J)**3
  620 CONTINUE
      DO 630 J = 1,5
          H(J* (J+1)/2) = SIGN(12.0D0,F1)*Y(85+J)*X(J)
  630 CONTINUE
      K = 1
      DO 650 I = 1,5
          DO 640 J = 1,I
              H(K) = H(K) + 2.0D0*Y(55+I+J*5)
              K = K + 1
  640     CONTINUE
  650 CONTINUE
      DO 670 J = 1,5
          T = -3.0D0*Y(85+J)*X(J)*X(J) - Y(90+J)
          DO 660 I = 1,15
              IF (I.LE.5) T = T - 2.0D0*Y(55+I+J*5)*X(I)
              IF (I.GT.5) T = T + Y(I+J*10-15)*X(I)
  660     CONTINUE
          IF (T.LE.0.0D0) GO TO 670
          I = J* (J+1)/2
          H(I) = H(I) - 6.0D2*Y(85+J)
  670 CONTINUE
      RETURN

      END
* SUBROUTINE TILD22             ALL SYSTEMS                 99/12/01
C PORTABILITY : ALL SYSTEMS
C 94/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
*  INITIATION OF VARIABLES FOR NONLINEAR MINIMAX APPROXIMATION.
*  LINEARLY CONSTRAINED DENSE VERSION.
*
* PARAMETERS :
*  IO  N  NUMBER OF VARIABLES.
*  IO  NA  NUMBER OF PARTIAL FUNCTIONS.
*  IO  NB  NUMBER OF BOX CONSTRAINTS.
*  IO  NC  NUMBER OF GENERAL LINEAR CONSTRAINTS.
*  RO  X(N)  VECTOR OF VARIABLES.
*  IO  IX(NF)  VECTOR CONTAINING TYPES OF BOUNDS.
*  RO  XL(NF)  VECTOR CONTAINING LOWER BOUNDS FOR VARIABLES.
*  RO  XU(NF)  VECTOR CONTAINING UPPER BOUNDS FOR VARIABLES.
*  IO  IC(NC)  VECTOR CONTAINING TYPES OF CONSTRAINTS.
*  RO  CL(NC)  VECTOR CONTAINING LOWER BOUNDS FOR CONSTRAINT FUNCTIONS.
*  RO  CU(NC)  VECTOR CONTAINING UPPER BOUNDS FOR CONSTRAINT FUNCTIONS.
*  RO  CG(NF*NC) MATRIX WHOSE COLUMNS ARE NORMALS OF THE LINEAR
*         CONSTRAINTS.
*  RO  FMIN  LOWER BOUND FOR VALUE OF THE OBJECTIVE FUNCTION.
*  RO  XMAX  MAXIMUM STEPSIZE.
*  IO  NEXT  NUMBER OF THE TEST PROBLEM.
*  IO  IEXT  TYPE OF OBJECTIVE FUNCTION. IEXT<0-MAXIMUM OF VALUES.
*         IEXT=0-MAXIMUM OF ABSOLUTE VALUES.
*  IO  IERR  ERROR INDICATOR.
*
      SUBROUTINE TILD22(N,NA,NB,NC,X,IX,XL,XU,IC,CL,CU,CG,FMIN,XMAX,
     +                  NEXT,IEXT,IERR)
C     .. Parameters ..
      DOUBLE PRECISION PI
      PARAMETER (PI=3.14159265358979323846D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION FMIN,XMAX
      INTEGER IERR,IEXT,N,NA,NB,NC,NEXT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CG(N*NC),CL(NC),CU(NC),X(N),XL(N),XU(N)
      INTEGER IC(NC),IX(N)
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION Y(163)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A
      INTEGER I,J,K,L
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC COS,DBLE,SIN
C     ..
C     .. Common blocks ..
      COMMON /EMPR22/Y
C     ..
      FMIN = -1.0D60
      XMAX = 1.0D3
      IEXT = -1
      IERR = 0
      NB = 0
      GO TO (10,20,30,40,50,100,160,200,220,
     +       240,260,310,340,380,420) NEXT

   10 IF (N.GE.2 .AND. NA.GE.3) THEN
          N = 2
          NA = 3
          NC = 1
          X(1) = 1.0D0
          X(2) = 2.0D0
          IC(1) = 1
          CL(1) = 0.5D0
          CG(1) = 1.0D0
          CG(2) = 1.0D0

      ELSE
          IERR = 1
      END IF

      RETURN

   20 IF (N.GE.2 .AND. NA.GE.3) THEN
          N = 2
          NA = 3
          NC = 1
          X(1) = -2.0D0
          X(2) = -1.0D0
          IC(1) = 2
          CU(1) = -2.5D0
          CG(1) = 3.0D0
          CG(2) = 1.0D0

      ELSE
          IERR = 1
      END IF

      RETURN

   30 IF (N.GE.2 .AND. NA.GE.3) THEN
          N = 2
          NA = 3
          NB = N
          NC = 1
          X(1) = -1.0D0
          X(2) = 1.0D-2
          IX(1) = 0
          IX(2) = 1
          XL(2) = 1.0D-2
          IC(1) = 1
          CL(1) = -5.0D-1
          CG(1) = 5.0D-2
          CG(2) = -1.0D0

      ELSE
          IERR = 1
      END IF

      RETURN

   40 IF (N.GE.2 .AND. NA.GE.3) THEN
          N = 2
          NA = 3
          NB = N
          NC = 1
          X(1) = -1.0D0
          X(2) = 3.0D0
          IX(1) = 0
          IX(2) = 1
          XL(2) = 1.0D-2
          IC(1) = 1
          CL(1) = 1.0D0
          CG(1) = -9.0D-1
          CG(2) = 1.0D0

      ELSE
          IERR = 1
      END IF

      RETURN

   50 IF (N.GE.6 .AND. NA.GE.3) THEN
          N = 6
          NA = 3
          X(1) = -1.0D0
          X(2) = 0.0D0
          X(3) = 0.0D0
          X(4) = -1.0D0
          X(5) = 1.0D0
          X(6) = 1.0D0
          NC = 5*NA
          DO 60 I = 1,NC
              CU(I) = 1.0D0
              IC(I) = 2
   60     CONTINUE
          DO 70 I = 1,N*NC
              CG(I) = 0.0D0
   70     CONTINUE
          K = 1
          DO 90 I = 1,NA
              L = 2* (I-1)
              DO 80 J = 1,5
                  CG(K+L) = SIN(2.0D0*PI*DBLE(J-1)/5.0D0)
                  CG(K+L+1) = COS(2.0D0*PI*DBLE(J-1)/5.0D0)
                  K = K + N
   80         CONTINUE
   90     CONTINUE

      ELSE
          IERR = 1
      END IF

      RETURN

  100 IF (N.GE.7 .AND. NA.GE.163) THEN
          N = 7
          NA = 163
          NB = N
          DO 110 I = 1,N
              X(I) = DBLE(I)*0.5D0
              IX(I) = 0
  110     CONTINUE
          XL(1) = 0.4D0
          IX(1) = 1
          IX(7) = 5
          DO 120 I = 1,NA
              Y(I) = 2.0D0*PI*SIN(PI* (8.5D0+DBLE(I)*0.5D0)/180.0D0)
  120     CONTINUE
          NC = 7
          DO 130 I = 1,6
              CL(I) = 0.4D0
              IC(I) = 1
  130     CONTINUE
          CL(7) = 1.0D0
          CU(7) = 1.0D0
          IC(7) = 5
          DO 140 I = 1,N*NC
              CG(I) = 0.0D0
  140     CONTINUE
          K = 0
          DO 150 I = 1,6
              CG(K+I) = -1.0D0
              CG(K+I+1) = 1.0D0
              K = K + N
  150     CONTINUE
          CG(46) = -1.0D0
          CG(48) = 1.0D0
          IEXT = 0
          FMIN = 0.0D0

      ELSE
          IERR = 1
      END IF

      RETURN

  160 IF (N.GE.8 .AND. NA.GE.8) THEN
          N = 8
          NA = 8
          NB = N
          DO 170 I = 1,N
              X(I) = 0.125D0
              XL(I) = 1.0D-8
              IX(I) = 1
  170     CONTINUE
          DO 180 I = 1,40
              Y(I) = 1.0D0
              Y(I+40) = 0.1D0
  180     CONTINUE
          Y(9) = 2.0D0
          Y(10) = 0.8D0
          Y(12) = 0.5D0
          Y(18) = 1.2D0
          Y(19) = 0.8D0
          Y(20) = 1.2D0
          Y(21) = 1.6D0
          Y(22) = 2.0D0
          Y(23) = 0.6D0
          Y(24) = 0.1D0
          Y(25) = 2.0D0
          Y(26) = 0.1D0
          Y(27) = 0.6D0
          Y(28) = 2.0D0
          Y(32) = 2.0D0
          Y(33) = 1.2D0
          Y(34) = 1.2D0
          Y(35) = 0.8D0
          Y(37) = 1.2D0
          Y(38) = 0.1D0
          Y(39) = 3.0D0
          Y(40) = 4.0D0
          Y(41) = 3.0D0
          Y(42) = 1.0D0
          Y(45) = 5.0D0
          Y(48) = 6.0D0
          Y(50) = 1.0D1
          Y(53) = 5.0D0
          Y(58) = 9.0D0
          Y(59) = 1.0D1
          Y(61) = 4.0D0
          Y(63) = 7.0D0
          Y(68) = 1.0D1
          Y(70) = 3.0D0
          Y(80) = 1.1D1
          Y(81) = 0.5D0
          Y(82) = 1.2D0
          Y(83) = 0.8D0
          Y(84) = 2.0D0
          Y(85) = 1.5D0
          NC = 1
          CL(1) = 1.0D0
          CU(1) = 1.0D0
          IC(1) = 5
          DO 190 I = 1,N
              CG(I) = 1.0D0
  190     CONTINUE

      ELSE
          IERR = 1
      END IF

      RETURN

  200 IF (N.GE.10 .AND. NA.GE.6) THEN
          N = 10
          NA = 6
          X(1) = 2.0D0
          X(2) = 3.0D0
          X(3) = 5.0D0
          X(4) = 5.0D0
          X(5) = 1.0D0
          X(6) = 2.0D0
          X(7) = 7.0D0
          X(8) = 3.0D0
          X(9) = 6.0D0
          X(10) = 1.0D1
          NC = 3
          CU(1) = 1.05D2
          CU(2) = 0.00D0
          CU(3) = 1.20D1
          IC(1) = 2
          IC(2) = 2
          IC(3) = 2
          DO 210 I = 1,N*NC
              CG(I) = 0.0D0
  210     CONTINUE
          CG(1) = 4.0D0
          CG(2) = 5.0D0
          CG(7) = -3.0D0
          CG(8) = 9.0D0
          CG(11) = 1.0D1
          CG(12) = -8.0D0
          CG(17) = -1.7D1
          CG(18) = 2.0D0
          CG(21) = -8.0D0
          CG(22) = 2.0D0
          CG(29) = 5.0D0
          CG(30) = -2.0D0

      ELSE
          IERR = 1
      END IF

      RETURN

  220 IF (N.GE.20 .AND. NA.GE.14) THEN
          N = 20
          NA = 14
          X(1) = 2.0D0
          X(2) = 3.0D0
          X(3) = 5.0D0
          X(4) = 5.0D0
          X(5) = 1.0D0
          X(6) = 2.0D0
          X(7) = 7.0D0
          X(8) = 3.0D0
          X(9) = 6.0D0
          X(10) = 1.0D1
          X(11) = 2.0D0
          X(12) = 2.0D0
          X(13) = 6.0D0
          X(14) = 1.5D1
          X(15) = 1.0D0
          X(16) = 2.0D0
          X(17) = 1.0D0
          X(18) = 2.0D0
          X(19) = 1.0D0
          X(20) = 3.0D0
          NC = 4
          CU(1) = 1.05D2
          CU(2) = 0.00D0
          CU(3) = 1.20D1
          CU(4) = 0.00D0
          IC(1) = 2
          IC(2) = 2
          IC(3) = 2
          IC(4) = 2
          DO 230 I = 1,N*NC
              CG(I) = 0.0D0
  230     CONTINUE
          CG(1) = 4.0D0
          CG(2) = 5.0D0
          CG(7) = -3.0D0
          CG(8) = 9.0D0
          CG(21) = 1.0D1
          CG(22) = -8.0D0
          CG(27) = -1.7D1
          CG(28) = 2.0D0
          CG(41) = -8.0D0
          CG(42) = 2.0D0
          CG(49) = 5.0D0
          CG(50) = -2.0D0
          CG(61) = 1.0D0
          CG(62) = 1.0D0
          CG(71) = 4.0D0
          CG(72) = -2.1D1

      ELSE
          IERR = 1
      END IF

      RETURN

  240 IF (N.GE.20 .AND. NA.GE.38) THEN
          N = 20
          NA = 38
          NB = N
          DO 250 I = 1,N
              X(I) = 1.0D2
              IX(I) = 0
              IF (I.LE.10) IX(I) = 1
              XL(I) = 0.5D0
  250     CONTINUE
          NC = 0
          IEXT = 0
          FMIN = 0.0D0

      ELSE
          IERR = 1
      END IF

      RETURN

  260 IF (N.GE.9 .AND. NA.GE.124) THEN
          N = 9
          NA = 124
          NB = N
          K = (N-1)/2
C      X(1)=1.8D-2
C      X(2)=1.9D-2
C      X(3)=2.0D-2
C      X(4)=2.1D-2
C      X(5)=0.8D 0
C      X(6)=0.9D 0
C      X(7)=1.0D 0
C      X(8)=1.1D 0
C      X(9)=-1.4D 1
          X(1) = 0.398D-1
          X(2) = 0.968D-4
          X(3) = 0.103D-3
          X(4) = 0.389D-1
          X(5) = 0.101D1
          X(6) = 0.968D0
          X(7) = 0.103D1
          X(8) = 0.988D0
          X(9) = -0.116D2
          DO 270 I = 1,N - 1
              XL(I) = 0.0D0
              IX(I) = 1
  270     CONTINUE
          IX(N) = 0
          L = (NA-2)/2
          A = (1.012577D0-0.987423D0)/DBLE(L-1)
          Y(1) = 0.967320D0
          DO 280 I = 2,L + 1
              Y(I) = 0.987423D0 + DBLE(I-2)*A
              Y(I+L) = Y(I)
  280     CONTINUE
          Y(NA) = 1.032680D0
          NC = K
          DO 290 I = 1,N*NC
              CG(I) = 0.0D0
  290     CONTINUE
          L = 0
          DO 300 I = 1,NC
              CG(L+I) = 1.0D4
              CG(L+I+K) = -1.0D0
              CL(I) = 0.0D0
              IC(I) = 1
              L = L + N
  300     CONTINUE
          IEXT = 0
          FMIN = 0.0D0

      ELSE
          IERR = 1
      END IF

      RETURN

  310 IF (N.GE.10 .AND. NA.GE.9) THEN
          N = 10
          NA = 9
          NB = N
          X(1) = 1745.0D0
          X(2) = 1200.0D1
          X(3) = 1100.0D-1
          X(4) = 3048.0D0
          X(5) = 1974.0D0
          X(6) = 8920.0D-2
          X(7) = 9280.0D-2
          X(8) = 8000.0D-3
          X(9) = 3600.0D-3
          X(10) = 145.0D0
          XL(1) = 1.0D-5
          XL(2) = 1.0D-5
          XL(3) = 1.0D-5
          XL(4) = 1.0D-5
          XL(5) = 1.0D-5
          XL(6) = 8.5D1
          XL(7) = 9.0D1
          XL(8) = 3.0D0
          XL(9) = 1.2D0
          XL(10) = 1.4D2
          XU(1) = 2.0D3
          XU(2) = 1.6D4
          XU(3) = 1.2D2
          XU(4) = 5.0D3
          XU(5) = 2.0D3
          XU(6) = 9.3D1
          XU(7) = 9.5D1
          XU(8) = 1.2D1
          XU(9) = 4.0D0
          XU(10) = 1.6D2
          DO 320 I = 1,N
              IX(I) = 3
  320     CONTINUE
          NC = 5
          CU(1) = 35.82D0
          CL(2) = 35.82D0
          CL(3) = 133.0D0
          CU(4) = 133.0D0
          CL(5) = 0.0D0
          CU(5) = 0.0D0
          IC(1) = 2
          IC(2) = 1
          IC(3) = 1
          IC(4) = 2
          IC(5) = 5
          DO 330 I = 1,N*NC
              CG(I) = 0.0D0
  330     CONTINUE
          CG(9) = 0.90D0
          CG(10) = 2.22D-1
          CG(19) = 1.00D0/0.90D0
          CG(20) = 2.22D-1
          CG(27) = 3.00D0
          CG(30) = -0.99D0
          CG(37) = 3.00D0
          CG(40) = -1.00D0/0.99D0
          CG(41) = -1.00D0
          CG(44) = 1.22D0
          CG(45) = -1.00D0

      ELSE
          IERR = 1
      END IF

      RETURN

  340 IF (N.GE.7 .AND. NA.GE.15) THEN
          N = 7
          NA = 13
          NB = N
          X(1) = 17.45D2
          X(2) = 1.10D2
          X(3) = 30.48D2
          X(4) = 8.90D1
          X(5) = 9.20D1
          X(6) = 8.00D0
          X(7) = 14.50D1
          XL(1) = 1.00D0
          XL(2) = 1.00D0
          XL(3) = 1.00D0
          XL(4) = 8.50D1
          XL(5) = 9.00D1
          XL(6) = 3.00D0
          XL(7) = 1.45D2
          XU(1) = 2.00D3
          XU(2) = 1.20D2
          XU(3) = 5.00D3
          XU(4) = 9.30D1
          XU(5) = 9.50D1
          XU(6) = 1.20D1
          XU(7) = 1.62D2
          DO 350 I = 1,N
              IX(I) = 3
  350     CONTINUE
          Y(1) = 1.71500000D0
          Y(2) = 0.03500000D0
          Y(3) = 4.05650000D0
          Y(4) = 10.0000000D0
          Y(5) = 3000.00000D0
          Y(6) = -0.06300000D0
          Y(7) = 0.59553571D-2
          Y(8) = 0.88392857D0
          Y(9) = -0.11756250D0
          Y(10) = 1.10880000D0
          Y(11) = 0.13035330D0
          Y(12) = -0.00660330D0
          Y(13) = 0.66173269D-3
          Y(14) = 0.17239878D-1
          Y(15) = -0.56595559D-2
          Y(16) = -0.19120592D-1
          Y(17) = 0.56850750D2
          Y(18) = 1.08702000D0
          Y(19) = 0.32175000D0
          Y(20) = -0.03762000D0
          Y(21) = 0.00619800D0
          Y(22) = 0.24623121D4
          Y(23) = -0.25125634D2
          Y(24) = 0.16118996D3
          Y(25) = 5000.00000D0
          Y(26) = -0.48951000D6
          Y(27) = 0.44333333D2
          Y(28) = 0.33000000D0
          Y(29) = 0.02255600D0
          Y(30) = -0.00759500D0
          Y(31) = 0.00061000D0
          Y(32) = -0.00050000D0
          Y(33) = 0.81967200D0
          Y(34) = 0.81967200D0
          Y(35) = 24500.0000D0
          Y(36) = -250.000000D0
          Y(37) = 0.10204082D-1
          Y(38) = 0.12244898D-4
          Y(39) = 0.00006250D0
          Y(40) = 0.00006250D0
          Y(41) = -0.00007625D0
          Y(42) = 1.22000000D0
          Y(43) = 1.00000000D0
          Y(44) = -1.00000000D0
          NC = 2
          L = 0
          DO 370 I = 1,NC
              CU(I) = 1.0D0
              IC(I) = 2
              DO 360 J = 1,N
                  CG(L+J) = 0.0D0
  360         CONTINUE
              L = L + N
  370     CONTINUE
          CG(5) = Y(29)
          CG(7) = Y(30)
          CG(8) = Y(32)
          CG(10) = Y(31)

      ELSE
          IERR = 1
      END IF

      RETURN

  380 IF (N.GE.8 .AND. NA.GE.7) THEN
          N = 8
          NA = 4
          NB = N
          X(1) = 5.00D3
          X(2) = 5.00D3
          X(3) = 5.00D3
          X(4) = 2.00D2
          X(5) = 3.50D2
          X(6) = 1.50D2
          X(7) = 2.25D2
          X(8) = 4.25D2
          DO 390 I = 1,N
              XL(I) = 1.0D1
              XU(I) = 1.0D3
              IX(I) = 3
  390     CONTINUE
          XL(1) = 1.0D2
          XL(2) = 1.0D3
          XL(3) = 1.0D3
          XU(1) = 1.0D4
          XU(2) = 1.0D4
          XU(3) = 1.0D4
          NC = 3
          L = 0
          DO 410 I = 1,NC
              CU(I) = 1.0D0
              IC(I) = 2
              DO 400 J = 1,N
                  CG(L+J) = 0.0D0
  400         CONTINUE
              L = L + N
  410     CONTINUE
          CG(4) = 2.5D-3
          CG(6) = 2.5D-3
          CG(12) = -2.5D-3
          CG(13) = 2.5D-3
          CG(15) = 2.5D-3
          CG(21) = -1.0D-2
          CG(24) = 1.0D-2

      ELSE
          IERR = 1
      END IF

      RETURN

  420 IF (N.GE.16 .AND. NA.GE.20) THEN
          N = 16
          NA = 19
          NB = N
          X(1) = 0.80D0
          X(2) = 0.83D0
          X(3) = 0.85D0
          X(4) = 0.87D0
          X(5) = 0.90D0
          X(6) = 0.10D0
          X(7) = 0.12D0
          X(8) = 0.19D0
          X(9) = 0.25D0
          X(10) = 0.29D0
          X(11) = 5.12D2
          X(12) = 1.31D1
          X(13) = 7.18D1
          X(14) = 6.40D2
          X(15) = 6.50D2
          X(16) = 5.70D0
          XL(1) = 1.0D-1
          XL(2) = 1.0D-1
          XL(3) = 1.0D-1
          XL(4) = 1.0D-1
          XL(5) = 9.0D-1
          XL(6) = 1.0D-4
          XL(7) = 1.0D-1
          XL(8) = 1.0D-1
          XL(9) = 1.0D-1
          XL(10) = 1.0D-1
          XL(11) = 1.0D0
          XL(12) = 1.0D-6
          XL(13) = 1.0D0
          XL(14) = 5.0D2
          XL(15) = 5.0D2
          XL(16) = 1.0D-6
          XU(1) = 9.0D-1
          XU(2) = 9.0D-1
          XU(3) = 9.0D-1
          XU(4) = 9.0D-1
          XU(5) = 1.0D0
          XU(6) = 1.0D-1
          XU(7) = 9.0D-1
          XU(8) = 9.0D-1
          XU(9) = 9.0D-1
          XU(10) = 9.0D-1
          XU(11) = 1.0D4
          XU(12) = 5.0D3
          XU(13) = 5.0D3
          XU(14) = 1.0D4
          XU(15) = 1.0D4
          XU(16) = 5.0D3
          DO 430 I = 1,N
              IX(I) = 3
  430     CONTINUE
          NC = 1
          CU(1) = 1.0D0
          IC(1) = 2
          DO 440 I = 1,N
              CG(I) = 0.0D0
  440     CONTINUE
          CG(11) = 2.0D-3
          CG(12) = -2.0D-3

      ELSE
          IERR = 1
      END IF

      RETURN

      END
* SUBROUTINE TAFU22             ALL SYSTEMS                 99/12/01
C PORTABILITY : ALL SYSTEMS
C 94/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
*  VALUES OF PARTIAL FUNCTIONS IN THE MINIMAX CRITERION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  II  KA  INDEX OF THE PARTIAL FUNCTION.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  FA  VALUE OF THE PARTIAL FUNCTION AT THE
*          SELECTED POINT.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*
      SUBROUTINE TAFU22(N,KA,X,FA,NEXT)
C     .. Scalar Arguments ..
      DOUBLE PRECISION FA
      INTEGER KA,N,NEXT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION X(N)
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION Y(163)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,B,P,S
      INTEGER I,J,K
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC COS,EXP,LOG,MOD,SIN,SINH,SQRT
C     ..
C     .. Common blocks ..
      COMMON /EMPR22/Y
C     ..
      GO TO (10,10,50,50,90,130,150,180,250,
     +       340,360,380,480,620,670) NEXT

   10 GO TO (20,30,40) KA

   20 FA = X(1)**2 + X(2)**2 + X(1)*X(2) - 1.0D0
      RETURN

   30 FA = SIN(X(1))
      RETURN

   40 FA = -COS(X(2))
      RETURN

   50 GO TO (60,70,80) KA

   60 FA = -EXP(X(1)-X(2))
      RETURN

   70 FA = SINH(X(1)-1.0D0) - 1.0D0
      RETURN

   80 FA = -LOG(X(2)) - 1.0D0
      RETURN

   90 GO TO (100,110,120) KA

  100 FA = -SQRT((X(1)-X(3))**2+ (X(2)-X(4))**2)
      RETURN

  110 FA = -SQRT((X(3)-X(5))**2+ (X(4)-X(6))**2)
      RETURN

  120 FA = -SQRT((X(5)-X(1))**2+ (X(6)-X(2))**2)
      RETURN

  130 A = 0.0D0
      DO 140 I = 1,N
          A = A + COS(Y(KA)*X(I))
  140 CONTINUE
      FA = (1.0D0+2.0D0*A)/1.5D1
      RETURN

  150 FA = 0.0D0
      K = 0
      DO 170 I = 1,5
          A = 0.0D0
          P = 0.0D0
          DO 160 J = 1,N
              A = A + Y(K+J)*X(J)** (1.0D0-Y(I+80))
              P = P + Y(K+J+40)*X(J)
  160     CONTINUE
          FA = FA + Y(K+KA)*P/ (X(KA)**Y(I+80)*A) - Y(K+KA+40)
          K = K + N
  170 CONTINUE
      RETURN

  180 FA = X(1)**2 + X(2)**2 + X(1)*X(2) - 1.4D1*X(1) - 1.6D1*X(2) +
     +     (X(3)-1.0D1)**2 + 4.0D0* (X(4)-5.0D0)**2 + (X(5)-3.0D0)**2 +
     +     2.0D0* (X(6)-1.0D0)**2 + 5.0D0*X(7)**2 +
     +     7.0D0* (X(8)-1.1D1)**2 + 2.0D0* (X(9)-1.0D1)**2 +
     +     (X(10)-7.0D0)**2 + 4.5D1
      GO TO (190,200,210,220,230,240) KA

  190 CONTINUE
      RETURN

  200 FA = FA + 1.0D1* (3.0D0* (X(1)-2.0D0)**2+4.0D0* (X(2)-3.0D0)**2+
     +     2.0D0*X(3)**2-7.0D0*X(4)-1.2D2)
      RETURN

  210 FA = FA + 1.0D1* (5.0D0*X(1)**2+8.0D0*X(2)+ (X(3)-6.0D0)**2-
     +     2.0D0*X(4)-4.0D1)
      RETURN

  220 FA = FA + 1.0D1* (0.5D0* (X(1)-8.0D0)**2+2.0D0* (X(2)-4.0D0)**2+
     +     3.0D0*X(5)**2-X(6)-3.0D1)
      RETURN

  230 FA = FA + 1.0D1* (X(1)**2+2.0D0* (X(2)-2.0D0)**2-2.0D0*X(1)*X(2)+
     +     1.4D1*X(5)-6.0D0*X(6))
      RETURN

  240 FA = FA + 1.0D1* (6.0D0*X(2)-3.0D0*X(1)+1.2D1* (X(9)-8.0D0)**2-
     +     7.0D0*X(10))
      RETURN

  250 FA = X(1)**2 + X(2)**2 + X(1)*X(2) - 1.4D1*X(1) - 1.6D1*X(2) +
     +     (X(3)-1.0D1)**2 + 4.0D0* (X(4)-5.0D0)**2 + (X(5)-3.0D0)**2 +
     +     2.0D0* (X(6)-1.0D0)**2 + 5.0D0*X(7)**2 +
     +     7.0D0* (X(8)-1.1D1)**2 + 2.0D0* (X(9)-1.0D1)**2 +
     +     (X(10)-7.0D0)**2 + (X(11)-9.0D0)**2 +
     +     1.0D1* (X(12)-1.0D0)**2 + 5.0D0* (X(13)-7.0D0)**2 +
     +     4.0D0* (X(14)-1.4D1)**2 + 2.7D1* (X(15)-1.0D0)**2 +
     +     X(16)**4 + (X(17)-2.0D0)**2 + 1.3D1* (X(18)-2.0D0)**2 +
     +     (X(19)-3.D0)**2 + X(20)**2 + 9.5D1
      GO TO (190,200,210,220,230,240,260,270,280,290,300,
     +       310,320,330) KA

  260 FA = FA + 1.0D1* (X(1)**2+1.5D1*X(11)-8.0D0*X(12)-2.8D1)
      RETURN

  270 FA = FA + 1.0D1* (4.0D0*X(1)+9.0D0*X(2)+5.0D0*X(13)**2-
     +     9.0D0*X(14)-8.7D1)
      RETURN

  280 FA = FA + 1.0D1* (3.0D0*X(1)+4.0D0*X(2)+3.0D0* (X(13)-6.0D0)**2-
     +     1.4D1*X(14)-1.0D1)
      RETURN

  290 FA = FA + 1.0D1* (1.4D1*X(1)**2+3.5D1*X(15)-7.9D1*X(16)-9.2D1)
      RETURN

  300 FA = FA + 1.0D1* (1.5D1*X(2)**2+1.1D1*X(15)-6.1D1*X(16)-5.4D1)
      RETURN

  310 FA = FA + 1.0D1* (5.0D0*X(1)**2+2.0D0*X(2)+9.0D0*X(17)**4-X(18)-
     +     6.8D1)
      RETURN

  320 FA = FA + 1.0D1* (X(1)**2-X(2)+1.9D1*X(19)-2.0D1*X(20)+1.9D1)
      RETURN

  330 FA = FA + 1.0D1* (7.0D0*X(1)**2+5.0D0*X(2)**2+X(19)**2-
     +     3.0D1*X(20))
      RETURN

  340 FA = -1.0D0
      DO 350 I = 1,N
          FA = FA + X(I)
  350 CONTINUE
      IF (MOD(KA,2).EQ.0) THEN
          I = (KA+2)/2
          FA = FA + X(I)* (X(I)-1.0D0)

      ELSE
          I = (KA+1)/2
          FA = FA + X(I)* (2.0D0*X(I)-1.0D0)
      END IF

      RETURN

  360 K = (N-1)/2
      A = Y(KA)
      S = 1.0D0
      IF (KA.GT.62 .AND. KA.LT.124) S = -S
      P = -8.0D0*LOG(A)
      A = A*A
      DO 370 I = 1,K
          B = X(I+K)**2 - A
          P = P + LOG(B*B+A*X(I)**2)
  370 CONTINUE
      FA = (0.5D0*P-X(N))*S
      IF (KA.EQ.1 .OR. KA.EQ.124) FA = FA + 3.0164D0
      RETURN

  380 P = 5.0D2
      A = 0.99D0
      FA = 5.04D0*X(1) + 0.35D-1*X(2) + 1.00D1*X(3) + 3.36D0*X(5) -
     +     0.63D-1*X(4)*X(7)
      GO TO (390,400,410,420,430,440,450,460,470) KA

  390 CONTINUE
      RETURN

  400 FA = FA + P* (1.12D0*X(1)+X(1)*X(8)* (1.3167D-1-6.67D-3*X(8))-
     +     X(4)/A)
      RETURN

  410 FA = FA - P* (1.12D0*X(1)+X(1)*X(8)* (1.3167D-1-6.67D-3*X(8))-
     +     X(4)*A)
      RETURN

  420 FA = FA + P* (57.425D0+X(8)* (1.098D0-0.038D0*X(8))+0.325D0*X(6)-
     +     X(7)/A)
      RETURN

  430 FA = FA - P* (57.425D0+X(8)* (1.098D0-0.038D0*X(8))+0.325D0*X(6)-
     +     X(7)*A)
      RETURN

  440 FA = FA + P* (9.8D4*X(3)/ (X(4)*X(9)+1.0D3*X(3))-X(6))
      RETURN

  450 FA = FA - P* (9.8D4*X(3)/ (X(4)*X(9)+1.0D3*X(3))-X(6))
      RETURN

  460 FA = FA + P* ((X(2)+X(5))/X(1)-X(8))
      RETURN

  470 FA = FA - P* ((X(2)+X(5))/X(1)-X(8))
      RETURN

  480 P = 1.0D5
      FA = Y(1)*X(1) + Y(2)*X(1)*X(6) + Y(3)*X(3) + Y(4)*X(2) + Y(5) +
     +     Y(6)*X(3)*X(5)
      GO TO (490,500,510,520,530,540,550,560,570,
     +       580,590,600,610) KA

  490 CONTINUE
      RETURN

  500 FA = FA + P* (Y(7)*X(6)**2+Y(8)*X(3)/X(1)+Y(9)*X(6)-1.0D0)
      RETURN

  510 FA = FA + P* ((Y(10)+Y(11)*X(6)+Y(12)*X(6)**2)*X(1)/X(3)-1.0D0)
      RETURN

  520 FA = FA + P* (Y(13)*X(6)**2+Y(14)*X(5)+Y(15)*X(4)+Y(16)*X(6)-
     +     1.0D0)
      RETURN

  530 FA = FA + P* ((Y(17)+Y(18)*X(6)+Y(19)*X(4)+Y(20)*X(6)**2)/X(5)-
     +     1.0D0)
      RETURN

  540 FA = FA + P* (Y(21)*X(7)+ (Y(22)/X(4)+Y(23))*X(2)/X(3)-1.0D0)
      RETURN

  550 FA = FA + P* ((Y(24)+ (Y(25)+Y(26)/X(4))*X(2)/X(3))/X(7)-1.0D0)
      RETURN

  560 FA = FA + P* ((Y(27)+Y(28)*X(7))/X(5)-1.0D0)
      RETURN

  570 FA = FA + P* ((Y(33)*X(1)+Y(34))/X(3)-1.0D0)
      RETURN

  580 FA = FA + P* ((Y(35)/X(4)+Y(36))*X(2)/X(3)-1.0D0)
      RETURN

  590 FA = FA + P* ((Y(37)+Y(38)*X(3)/X(2))*X(4)-1.0D0)
      RETURN

  600 FA = FA + P* (Y(39)*X(1)*X(6)+Y(40)*X(1)+Y(41)*X(3)-1.0D0)
      RETURN

  610 FA = FA + P* ((Y(42)*X(3)+Y(43))/X(1)+Y(44)*X(6)-1.0D0)
      RETURN

  620 P = 1.0D5
      FA = X(1) + X(2) + X(3)
      GO TO (630,640,650,660) KA

  630 CONTINUE
      RETURN

  640 FA = FA + P* ((833.33252D0*X(4)/X(1)+1.0D2-83333.333D0/X(1))/X(6)-
     +     1.0D0)
      RETURN

  650 FA = FA + P* ((1.25D3* (X(5)-X(4))/X(2)+X(4))/X(7)-1.0D0)
      RETURN

  660 FA = FA + P* (((1.25D6-2.5D3*X(5))/X(3)+X(5))/X(8)-1.0D0)
      RETURN

  670 P = 2.0D3
      FA = 1.262626D0* (X(12)+X(13)+X(14)+X(15)+X(16)) -
     +     1.231060D0* (X(1)*X(12)+X(2)*X(13)+X(3)*X(14)+X(4)*X(15)+
     +     X(5)*X(16))
      GO TO (680,690,700,710,720,730,740,750,760,
     +       770,780,790,800,810,820,830,840,850,
     +       860) KA

  680 CONTINUE
      RETURN

  690 FA = FA + P* (X(1)* (9.75D-1+ (3.475D-2-9.75D-3*X(1))/X(6))-1.0D0)
      RETURN

  700 FA = FA + P* (X(2)* (9.75D-1+ (3.475D-2-9.75D-3*X(2))/X(7))-1.0D0)
      RETURN

  710 FA = FA + P* (X(3)* (9.75D-1+ (3.475D-2-9.75D-3*X(3))/X(8))-1.0D0)
      RETURN

  720 FA = FA + P* (X(4)* (9.75D-1+ (3.475D-2-9.75D-3*X(4))/X(9))-1.0D0)
      RETURN

  730 FA = FA + P* (X(5)* (9.75D-1+ (3.475D-2-9.75D-3*X(5))/X(10))-
     +     1.0D0)
      RETURN

  740 FA = FA + P* ((X(6)+ (X(1)-X(6))*X(12)/X(11))/X(7)-1.0D0)
      RETURN

  750 FA = FA + P* (((X(7)+2.0D-3* (X(7)-X(1))*X(12))/X(8)+
     +     2.0D-3*X(13)* (X(2)/X(8)-1.0D0))-1.0D0)
      RETURN

  760 FA = FA + P* (X(8)+X(9)+2.0D-3*X(13)* (X(8)-X(2))+
     +     2.0D-3*X(14)* (X(3)-X(9))-1.0D0)
      RETURN

  770 FA = FA + P* ((X(9)+ ((X(4)-X(8))*X(15)+5.0D2* (X(10)-
     +     X(9)))/X(14))/X(3)-1.0D0)
      RETURN

  780 FA = FA + P* ((X(10)/X(4)+ (X(5)/X(4)-1.0D0)*X(16)/X(15)+
     +     5.0D2* (1.0D0-X(10)/X(4))/X(15))-1.0D0)
      RETURN

  790 FA = FA + P* (9.0D-1/X(4)+2.0D-3*X(16)* (1.0D0-X(5)/X(4))-1.0D0)
      RETURN

  800 FA = FA + P* (X(12)/X(11)-1.0D0)
      RETURN

  810 FA = FA + P* (X(4)/X(5)-1.0D0)
      RETURN

  820 FA = FA + P* (X(3)/X(4)-1.0D0)
      RETURN

  830 FA = FA + P* (X(2)/X(3)-1.0D0)
      RETURN

  840 FA = FA + P* (X(1)/X(2)-1.0D0)
      RETURN

  850 FA = FA + P* (X(9)/X(10)-1.0D0)
      RETURN

  860 FA = FA + P*X(8)/X(9) - P
      RETURN

      END
* SUBROUTINE TAGU22             ALL SYSTEMS                 99/12/01
C PORTABILITY : ALL SYSTEMS
C 94/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
*  GRADIENTS OF PARTIAL FUNCTIONS IN THE MINIMAX CRITERION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  II  KA  INDEX OF THE PARTIAL FUNCTION.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  GA(N)  GRADIENT OF THE PARTIAL FUNCTION AT THE
*          SELECTED POINT.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*
      SUBROUTINE TAGU22(N,KA,X,GA,NEXT)
C     .. Scalar Arguments ..
      INTEGER KA,N,NEXT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION GA(N),X(N)
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION Y(163)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,B,C,P,S
      INTEGER I,J,K
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC COS,COSH,EXP,MOD,SIN,SQRT
C     ..
C     .. Common blocks ..
      COMMON /EMPR22/Y
C     ..
      GO TO (10,10,50,50,90,130,150,200,270,
     +       360,380,400,500,640,690) NEXT

   10 GO TO (20,30,40) KA

   20 GA(1) = 2.0D0*X(1) + X(2)
      GA(2) = 2.0D0*X(2) + X(1)
      RETURN

   30 GA(1) = COS(X(1))
      GA(2) = 0.0D0
      RETURN

   40 GA(1) = 0.0D0
      GA(2) = SIN(X(2))
      RETURN

   50 GO TO (60,70,80) KA

   60 GA(1) = -EXP(X(1)-X(2))
      GA(2) = EXP(X(1)-X(2))
      RETURN

   70 GA(1) = COSH(X(1)-1.0D0)
      GA(2) = 0.0D0
      RETURN

   80 GA(1) = 0.0D0
      GA(2) = -1.0D0/X(2)
      RETURN

   90 GO TO (100,110,120) KA

  100 A = SQRT((X(1)-X(3))**2+ (X(2)-X(4))**2)
      GA(1) = - (X(1)-X(3))/A
      GA(2) = - (X(2)-X(4))/A
      GA(3) = -GA(1)
      GA(4) = -GA(2)
      GA(5) = 0.0D0
      GA(6) = 0.0D0
      RETURN

  110 A = SQRT((X(3)-X(5))**2+ (X(4)-X(6))**2)
      GA(1) = 0.0D0
      GA(2) = 0.0D0
      GA(3) = - (X(3)-X(5))/A
      GA(4) = - (X(4)-X(6))/A
      GA(5) = -GA(3)
      GA(6) = -GA(4)
      RETURN

  120 A = SQRT((X(5)-X(1))**2+ (X(6)-X(2))**2)
      GA(1) = (X(5)-X(1))/A
      GA(2) = (X(6)-X(2))/A
      GA(3) = 0.0D0
      GA(4) = 0.0D0
      GA(5) = -GA(1)
      GA(6) = -GA(2)
      RETURN

  130 DO 140 I = 1,N
          GA(I) = -2.0D0*Y(KA)*SIN(Y(KA)*X(I))/1.5D1
  140 CONTINUE
      RETURN

  150 DO 160 I = 1,N
          GA(I) = 0.0D0
  160 CONTINUE
      K = 0
      DO 190 I = 1,5
          A = 0.0D0
          P = 0.0D0
          DO 170 J = 1,N
              A = A + Y(K+J)*X(J)** (1.0D0-Y(I+80))
              P = P + Y(K+J+40)*X(J)
  170     CONTINUE
          B = Y(K+KA)/ (X(KA)**Y(I+80)*A)
          DO 180 J = 1,N
              C = Y(K+J)* (1.0D0-Y(I+80))/ (X(J)**Y(I+80)*A)
              GA(J) = GA(J) + B* (Y(K+J+40)-C*P)
  180     CONTINUE
          GA(KA) = GA(KA) - B*Y(I+80)*P/X(KA)
          K = K + N
  190 CONTINUE
      RETURN

  200 GA(1) = 2.0D0*X(1) + X(2) - 1.4D1
      GA(2) = 2.0D0*X(2) + X(1) - 1.6D1
      GA(3) = 2.0D0* (X(3)-1.0D1)
      GA(4) = 8.0D0* (X(4)-5.0D0)
      GA(5) = 2.0D0* (X(5)-3.0D0)
      GA(6) = 4.0D0* (X(6)-1.0D0)
      GA(7) = 1.0D1*X(7)
      GA(8) = 1.4D1* (X(8)-1.1D1)
      GA(9) = 4.0D0* (X(9)-1.0D1)
      GA(10) = 2.0D0* (X(10)-7.0D0)
      GO TO (210,220,230,240,250,260) KA

  210 CONTINUE
      RETURN

  220 GA(1) = GA(1) + 6.0D1* (X(1)-2.0D0)
      GA(2) = GA(2) + 8.0D1* (X(2)-3.0D0)
      GA(3) = GA(3) + 4.0D1*X(3)
      GA(4) = GA(4) - 7.0D1
      RETURN

  230 GA(1) = GA(1) + 1.0D2*X(1)
      GA(2) = GA(2) + 8.0D1
      GA(3) = GA(3) + 2.0D1* (X(3)-6.0D0)
      GA(4) = GA(4) - 2.0D1
      RETURN

  240 GA(1) = GA(1) + 1.0D1* (X(1)-8.0D0)
      GA(2) = GA(2) + 4.0D1* (X(2)-4.0D0)
      GA(5) = GA(5) + 6.0D1*X(5)
      GA(6) = GA(6) - 1.0D1
      RETURN

  250 GA(1) = GA(1) + 2.0D1*X(1) - 2.0D1*X(2)
      GA(2) = GA(2) + 4.0D1* (X(2)-2.0D0) - 2.0D1*X(1)
      GA(5) = GA(5) + 1.4D2
      GA(6) = GA(6) - 6.0D1
      RETURN

  260 GA(1) = GA(1) - 3.0D1
      GA(2) = GA(2) + 6.0D1
      GA(9) = GA(9) + 2.4D2* (X(9)-8.0D0)
      GA(10) = GA(10) - 7.0D1
      RETURN

  270 GA(1) = 2.0D0*X(1) + X(2) - 1.4D1
      GA(2) = 2.0D0*X(2) + X(1) - 1.6D1
      GA(3) = 2.0D0* (X(3)-1.0D1)
      GA(4) = 8.0D0* (X(4)-5.0D0)
      GA(5) = 2.0D0* (X(5)-3.0D0)
      GA(6) = 4.0D0* (X(6)-1.0D0)
      GA(7) = 1.0D1*X(7)
      GA(8) = 1.4D1* (X(8)-1.1D1)
      GA(9) = 4.0D0* (X(9)-1.0D1)
      GA(10) = 2.0D0* (X(10)-7.0D0)
      GA(11) = 2.0D0* (X(11)-9.0D0)
      GA(12) = 2.0D1* (X(12)-1.0D0)
      GA(13) = 1.0D1* (X(13)-7.0D0)
      GA(14) = 8.0D0* (X(14)-1.4D1)
      GA(15) = 5.4D1* (X(15)-1.0D0)
      GA(16) = 4.0D0*X(16)**3
      GA(17) = 2.0D0* (X(17)-2.0D0)
      GA(18) = 2.6D1* (X(18)-2.0D0)
      GA(19) = 2.0D0* (X(19)-3.0D0)
      GA(20) = 2.0D0*X(20)
      GO TO (210,220,230,240,250,260,280,290,300,310,320,
     +       330,340,350) KA

  280 GA(1) = GA(1) + 2.0D1*X(1)
      GA(11) = GA(11) + 1.5D2
      GA(12) = GA(12) - 8.0D1
      RETURN

  290 GA(1) = GA(1) + 4.0D1
      GA(2) = GA(2) + 9.0D1
      GA(13) = GA(13) + 1.0D2*X(13)
      GA(14) = GA(14) - 9.0D1
      RETURN

  300 GA(1) = GA(1) + 3.0D1
      GA(2) = GA(2) + 4.0D1
      GA(13) = GA(13) + 6.0D1* (X(13)-6.0D0)
      GA(14) = GA(14) - 1.4D2
      RETURN

  310 GA(1) = GA(1) + 2.8D2*X(1)
      GA(15) = GA(15) + 3.5D2
      GA(16) = GA(16) - 7.9D2
      RETURN

  320 GA(2) = GA(2) + 3.0D2*X(2)
      GA(15) = GA(15) + 1.1D2
      GA(16) = GA(16) - 6.1D2
      RETURN

  330 GA(1) = GA(1) + 1.0D2*X(1)
      GA(2) = GA(2) + 2.0D1
      GA(17) = GA(17) + 3.6D2*X(17)**3
      GA(18) = GA(18) - 1.0D1
      RETURN

  340 GA(1) = GA(1) + 2.0D1*X(1)
      GA(2) = GA(2) - 1.0D1
      GA(19) = GA(19) + 1.9D2
      GA(20) = GA(20) - 2.0D2
      RETURN

  350 GA(1) = GA(1) + 1.4D2*X(1)
      GA(2) = GA(2) + 1.0D2*X(2)
      GA(19) = GA(19) + 2.0D1*X(19)
      GA(20) = GA(20) - 3.0D2
      RETURN

  360 DO 370 I = 1,N
          GA(I) = 1.0D0
  370 CONTINUE
      IF (MOD(KA,2).EQ.0) THEN
          I = (KA+2)/2
          GA(I) = GA(I) + 2.0D0*X(I) - 1.0D0

      ELSE
          I = (KA+1)/2
          GA(I) = GA(I) + 4.0D0*X(I) - 1.0D0
      END IF

      RETURN

  380 K = (N-1)/2
      A = Y(KA)**2
      S = 1.0D0
      IF (KA.GT.62 .AND. KA.LT.124) S = -S
      DO 390 I = 1,K
          B = X(I+K)**2 - A
          P = S* (B*B+A*X(I)**2)
          GA(I) = A*X(I)/P
          GA(I+K) = 2.0D0*X(I+K)*B/P
  390 CONTINUE
      GA(N) = -S
      RETURN

  400 P = 5.0D2
      A = 0.99D0
      GA(1) = 5.04D0
      GA(2) = 0.35D-1
      GA(3) = 1.00D1
      GA(4) = -0.63D-1*X(7)
      GA(5) = 3.36D0
      GA(6) = 0.00D0
      GA(7) = -0.63D-1*X(4)
      GA(8) = 0.00D0
      GA(9) = 0.00D0
      GA(10) = 0.00D0
      GO TO (410,420,430,440,450,460,470,480,490) KA

  410 CONTINUE
      RETURN

  420 GA(1) = GA(1) + P* (1.12D0+X(8)* (1.3167D-1-6.67D-3*X(8)))
      GA(4) = GA(4) - P/A
      GA(8) = GA(8) + P*X(1)* (1.3167D-1-1.334D-2*X(8))
      RETURN

  430 GA(1) = GA(1) - P* (1.12D0+X(8)* (1.3167D-1-6.67D-3*X(8)))
      GA(4) = GA(4) + P*A
      GA(8) = GA(8) - P*X(1)* (1.3167D-1-1.334D-2*X(8))
      RETURN

  440 GA(6) = GA(6) + P*0.325D0
      GA(7) = GA(7) - P/A
      GA(8) = GA(8) + P* (1.098D0-0.076D0*X(8))
      RETURN

  450 GA(6) = GA(6) - P*0.325D0
      GA(7) = GA(7) + P*A
      GA(8) = GA(8) - P* (1.098D0-0.076D0*X(8))
      RETURN

  460 C = (X(4)*X(9)+1.0D3*X(3))**2
      GA(3) = GA(3) + 9.8D4*P*X(4)*X(9)/C
      GA(4) = GA(4) - 9.8D4*P*X(3)*X(9)/C
      GA(6) = GA(6) - P
      GA(9) = GA(9) - 9.8D4*P*X(3)*X(4)/C
      RETURN

  470 C = (X(4)*X(9)+1.0D3*X(3))**2
      GA(3) = GA(3) - 9.8D4*P*X(4)*X(9)/C
      GA(4) = GA(4) + 9.8D4*P*X(3)*X(9)/C
      GA(6) = GA(6) + P
      GA(9) = GA(9) + 9.8D4*P*X(3)*X(4)/C
      RETURN

  480 GA(1) = GA(1) - P* (X(2)+X(5))/X(1)**2
      GA(2) = GA(2) + P/X(1)
      GA(5) = GA(5) + P/X(1)
      GA(8) = GA(8) - P
      RETURN

  490 GA(1) = GA(1) + P* (X(2)+X(5))/X(1)**2
      GA(2) = GA(2) - P/X(1)
      GA(5) = GA(5) - P/X(1)
      GA(8) = GA(8) + P
      RETURN

  500 P = 1.0D5
      GA(1) = Y(1) + Y(2)*X(6)
      GA(2) = Y(4)
      GA(3) = Y(3) + Y(6)*X(5)
      GA(4) = 0.0D0
      GA(5) = Y(6)*X(3)
      GA(6) = Y(2)*X(1)
      GA(7) = 0.0D0
      GO TO (510,520,530,540,550,560,570,580,590,
     +       600,610,620,630) KA

  510 CONTINUE
      RETURN

  520 GA(1) = GA(1) - P*Y(8)*X(3)/X(1)**2
      GA(3) = GA(3) + P*Y(8)/X(1)
      GA(6) = GA(6) + P* (2.0D0*Y(7)*X(6)+Y(9))
      RETURN

  530 GA(1) = GA(1) + P* (Y(10)+Y(11)*X(6)+Y(12)*X(6)**2)/X(3)
      GA(3) = GA(3) - P* (Y(10)+Y(11)*X(6)+Y(12)*X(6)**2)*X(1)/X(3)**2
      GA(6) = GA(6) + P* (Y(11)+2.0D0*Y(12)*X(6))*X(1)/X(3)
      RETURN

  540 GA(4) = GA(4) + P*Y(15)
      GA(5) = GA(5) + P*Y(14)
      GA(6) = GA(6) + P* (2.0D0*Y(13)*X(6)+Y(16))
      RETURN

  550 GA(4) = GA(4) + P*Y(19)/X(5)
      GA(5) = GA(5) - P* (Y(17)+Y(18)*X(6)+Y(19)*X(4)+Y(20)*X(6)**2)/
     +        X(5)**2
      GA(6) = GA(6) + P* (Y(18)+2.0D0*Y(20)*X(6))/X(5)
      RETURN

  560 GA(2) = GA(2) + P* (Y(22)/X(4)+Y(23))/X(3)
      GA(3) = GA(3) - P* (Y(22)/X(4)+Y(23))*X(2)/X(3)**2
      GA(4) = GA(4) - P*Y(22)*X(2)/ (X(3)*X(4)**2)
      GA(7) = GA(7) + P*Y(21)
      RETURN

  570 GA(2) = GA(2) + P* (Y(25)+Y(26)/X(4))/ (X(3)*X(7))
      GA(3) = GA(3) - P* (Y(25)+Y(26)/X(4))*X(2)/ (X(3)**2*X(7))
      GA(4) = GA(4) - P*Y(26)*X(2)/ (X(3)*X(4)**2*X(7))
      GA(7) = GA(7) - P* (Y(24)+ (Y(25)+Y(26)/X(4))*X(2)/X(3))/X(7)**2
      RETURN

  580 GA(5) = GA(5) - P* (Y(27)+Y(28)*X(7))/X(5)**2
      GA(7) = GA(7) + P*Y(28)/X(5)
      RETURN

  590 GA(1) = GA(1) + P*Y(33)/X(3)
      GA(3) = GA(3) - P* (Y(33)*X(1)+Y(34))/X(3)**2
      RETURN

  600 GA(2) = GA(2) + P* (Y(35)/X(4)+Y(36))/X(3)
      GA(3) = GA(3) - P* (Y(35)/X(4)+Y(36))*X(2)/X(3)**2
      GA(4) = GA(4) - P*Y(35)*X(2)/ (X(3)*X(4)**2)
      RETURN

  610 GA(2) = GA(2) - P*Y(38)*X(3)*X(4)/X(2)**2
      GA(3) = GA(3) + P*Y(38)*X(4)/X(2)
      GA(4) = GA(4) + P* (Y(37)+Y(38)*X(3)/X(2))
      RETURN

  620 GA(1) = GA(1) + P* (Y(39)*X(6)+Y(40))
      GA(3) = GA(3) + P*Y(41)
      GA(6) = GA(6) + P*Y(39)*X(1)
      RETURN

  630 GA(1) = GA(1) - P* (Y(42)*X(3)+Y(43))/X(1)**2
      GA(3) = GA(3) + P*Y(42)/X(1)
      GA(6) = GA(6) + P*Y(44)
      RETURN

  640 P = 1.0D5
      GA(1) = 1.0D0
      GA(2) = 1.0D0
      GA(3) = 1.0D0
      GA(4) = 0.0D0
      GA(5) = 0.0D0
      GA(6) = 0.0D0
      GA(7) = 0.0D0
      GA(8) = 0.0D0
      GO TO (650,660,670,680) KA

  650 CONTINUE
      RETURN

  660 GA(1) = GA(1) - P* (833.33252D0*X(4)-83333.333D0)/ (X(1)**2*X(6))
      GA(4) = GA(4) + P*833.33252D0/ (X(1)*X(6))
      GA(6) = GA(6) - P* (833.33252D0*X(4)/X(1)+1.0D2-83333.333D0/X(1))/
     +        X(6)**2
      RETURN

  670 GA(2) = GA(2) - P*1.25D3* (X(5)-X(4))/ (X(2)**2*X(7))
      GA(4) = GA(4) + P* (1.0D0-1.25D3/X(2))/X(7)
      GA(5) = GA(5) + P*1.25D3/ (X(2)*X(7))
      GA(7) = GA(7) - P* (1.25D3* (X(5)-X(4))/X(2)+X(4))/X(7)**2
      RETURN

  680 GA(3) = GA(3) - P* (1.25D6-2.5D3*X(5))/ (X(3)**2*X(8))
      GA(5) = GA(5) + P* (1.0D0-2.5D3/X(3))/X(8)
      GA(8) = GA(8) - P* ((1.25D6-2.5D3*X(5))/X(3)+X(5))/X(8)**2
      RETURN

  690 P = 2.0D3
      GA(1) = -1.231060D0*X(12)
      GA(2) = -1.231060D0*X(13)
      GA(3) = -1.231060D0*X(14)
      GA(4) = -1.231060D0*X(15)
      GA(5) = -1.231060D0*X(16)
      GA(6) = 0.0D0
      GA(7) = 0.0D0
      GA(8) = 0.0D0
      GA(9) = 0.0D0
      GA(10) = 0.0D0
      GA(11) = 0.0D0
      GA(12) = 1.262626D0 - 1.231060D0*X(1)
      GA(13) = 1.262626D0 - 1.231060D0*X(2)
      GA(14) = 1.262626D0 - 1.231060D0*X(3)
      GA(15) = 1.262626D0 - 1.231060D0*X(4)
      GA(16) = 1.262626D0 - 1.231060D0*X(5)
      GO TO (700,710,720,730,740,750,760,770,780,
     +       790,800,810,820,830,840,850,860,870,
     +       880) KA

  700 CONTINUE
      RETURN

  710 GA(1) = GA(1) + P* (9.75D-1+ (3.475D-2-1.95D-2*X(1))/X(6))
      GA(6) = GA(6) - P*X(1)* (3.475D-2-9.75D-3*X(1))/X(6)**2
      RETURN

  720 GA(2) = GA(2) + P* (9.75D-1+ (3.475D-2-1.95D-2*X(2))/X(7))
      GA(7) = GA(7) - P*X(2)* (3.475D-2-9.75D-3*X(2))/X(7)**2
      RETURN

  730 GA(3) = GA(3) + P* (9.75D-1+ (3.475D-2-1.95D-2*X(3))/X(8))
      GA(8) = GA(8) - P*X(3)* (3.475D-2-9.75D-3*X(3))/X(8)**2
      RETURN

  740 GA(4) = GA(4) + P* (9.75D-1+ (3.475D-2-1.95D-2*X(4))/X(9))
      GA(9) = GA(9) - P*X(4)* (3.475D-2-9.75D-3*X(4))/X(9)**2
      RETURN

  750 GA(5) = GA(5) + P* (9.75D-1+ (3.475D-2-1.95D-2*X(5))/X(10))
      GA(10) = GA(10) - P*X(5)* (3.475D-2-9.75D-3*X(5))/X(10)**2
      RETURN

  760 GA(1) = GA(1) + P*X(12)/ (X(7)*X(11))
      GA(6) = GA(6) + P* (1.0D0-X(12)/X(11))/X(7)
      GA(7) = GA(7) - P* (X(6)+ (X(1)-X(6))*X(12)/X(11))/X(7)**2
      GA(11) = GA(11) - P* (X(1)-X(6))*X(12)/ (X(7)*X(11)**2)
      GA(12) = GA(12) + P* (X(1)-X(6))/ (X(11)*X(7))
      RETURN

  770 GA(1) = GA(1) - P*2.0D-3*X(12)/X(8)
      GA(2) = GA(2) + P*2.0D-3*X(13)/X(8)
      GA(7) = GA(7) + P* (1.0D0+2.0D-3*X(12))/X(8)
      GA(8) = GA(8) - P* (X(7)+2.0D-3* ((X(7)-X(1))*X(12)+X(2)*X(13)))/
     +        X(8)**2
      GA(12) = GA(12) + P*2.0D-3* (X(7)-X(1))/X(8)
      GA(13) = GA(13) + P*2.0D-3* (X(2)/X(8)-1.0D0)
      RETURN

  780 GA(2) = GA(2) - P*2.0D-3*X(13)
      GA(3) = GA(3) + P*2.0D-3*X(14)
      GA(8) = GA(8) + P* (1.0D0+2.0D-3*X(13))
      GA(9) = GA(9) + P* (1.0D0-2.0D-3*X(14))
      GA(13) = GA(13) + P*2.0D-3* (X(8)-X(2))
      GA(14) = GA(14) + P*2.0D-3* (X(3)-X(9))
      RETURN

  790 GA(3) = GA(3) - P* (X(9)+ ((X(4)-X(8))*X(15)+5.0D2* (X(10)-X(9)))/
     +        X(14))/X(3)**2
      GA(4) = GA(4) + P*X(15)/ (X(3)*X(14))
      GA(8) = GA(8) - P*X(15)/ (X(3)*X(14))
      GA(9) = GA(9) + P* (1.0D0-5.0D2/X(14))/X(3)
      GA(10) = GA(10) + P*5.0D2/ (X(3)*X(14))
      GA(14) = GA(14) - P* ((X(4)-X(8))*X(15)+5.0D2* (X(10)-X(9)))/
     +         (X(3)*X(14)**2)
      GA(15) = GA(15) + P* (X(4)-X(8))/ (X(3)*X(14))
      RETURN

  800 GA(4) = GA(4) - P* (X(10)+ (X(5)*X(16)-5.0D2*X(10))/X(15))/X(4)**2
      GA(5) = GA(5) + P*X(16)/ (X(4)*X(15))
      GA(10) = GA(10) + P* (1.0D0-5.0D2/X(15))/X(4)
      GA(15) = GA(15) - P* (5.0D2-X(16)+ (X(5)*X(16)-5.0D2*X(10))/X(4))/
     +         X(15)**2
      GA(16) = GA(16) + P* (X(5)/X(4)-1.0D0)/X(15)
      RETURN

  810 GA(4) = GA(4) - P* (9.0D-1-2.0D-3*X(5)*X(16))/X(4)**2
      GA(5) = GA(5) - P*2.0D-3*X(16)/X(4)
      GA(16) = GA(16) + P*2.0D-3* (1.0D0-X(5)/X(4))
      RETURN

  820 GA(11) = GA(11) - P*X(12)/X(11)**2
      GA(12) = GA(12) + P/X(11)
      RETURN

  830 GA(4) = GA(4) + P/X(5)
      GA(5) = GA(5) - P*X(4)/X(5)**2
      RETURN

  840 GA(3) = GA(3) + P/X(4)
      GA(4) = GA(4) - P*X(3)/X(4)**2
      RETURN

  850 GA(2) = GA(2) + P/X(3)
      GA(3) = GA(3) - P*X(2)/X(3)**2
      RETURN

  860 GA(1) = GA(1) + P/X(2)
      GA(2) = GA(2) - P*X(1)/X(2)**2
      RETURN

  870 GA(9) = GA(9) + P/X(10)
      GA(10) = GA(10) - P*X(9)/X(10)**2
      RETURN

  880 GA(8) = GA(8) + P/X(9)
      GA(9) = GA(9) - P*X(8)/X(9)**2
      RETURN

      END
* SUBROUTINE TAHD22             ALL SYSTEMS                 99/12/01
C PORTABILITY : ALL SYSTEMS
C 95/12/01 LU : ORIGINAL VERSION
*
* PURPOSE :
*  HESSIAN MATRICES OF PARTIAL FUNCTIONS IN THE MINIMAX CRITERION.
*  DENSE VERSION.
*
* PARAMETERS :
*  II  N  NUMBER OF VARIABLES.
*  II  KA  INDEX OF THE PARTIAL FUNCTION.
*  RI  X(N)  VECTOR OF VARIABLES.
*  RO  HA(N*(N+1)/2)  GRADIENT OF THE PARTIAL FUNCTION
*         AT THE SELECTED POINT.
*  II  NEXT  NUMBER OF THE TEST PROBLEM.
*
      SUBROUTINE TAHD22(N,KA,X,HA,NEXT)
C     .. Scalar Arguments ..
      INTEGER KA,N,NEXT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION HA(N* (N+1)/2),X(N)
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION Y(163)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,B,C,P,Q,R,S
      INTEGER I,J,K,KK,L,LL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC COS,EXP,MOD,SIN,SINH,SQRT
C     ..
C     .. Common blocks ..
      COMMON /EMPR22/Y
C     ..
      GO TO (10,10,50,50,90,140,170,230,310,
     +       410,430,460,570,710,770) NEXT

   10 GO TO (20,30,40) KA

   20 HA(1) = 2.0D0
      HA(2) = 1.0D0
      HA(3) = 2.0D0
      RETURN

   30 HA(1) = -SIN(X(1))
      HA(2) = 0.0D0
      HA(3) = 0.0D0
      RETURN

   40 HA(1) = 0.0D0
      HA(2) = 0.0D0
      HA(3) = COS(X(2))
      RETURN

   50 GO TO (60,70,80) KA

   60 HA(1) = -EXP(X(1)-X(2))
      HA(2) = EXP(X(1)-X(2))
      HA(3) = -EXP(X(1)-X(2))
      RETURN

   70 HA(1) = SINH(X(1)-1.0D0)
      HA(2) = 0.0D0
      HA(3) = 0.0D0
      RETURN

   80 HA(1) = 0.0D0
      HA(2) = 0.0D0
      HA(3) = 1.0D0/X(2)**2
      RETURN

   90 DO 100 I = 1,N* (N+1)/2
          HA(I) = 0.0D0
  100 CONTINUE
      GO TO (110,120,130) KA

  110 A = SQRT((X(1)-X(3))**2+ (X(2)-X(4))**2)
      B = (X(1)-X(3))/A
      C = (X(2)-X(4))/A
      HA(1) = -1.0D0/A + B*B/A
      HA(2) = B*C/A
      HA(3) = -1.0D0/A + C*C/A
      HA(4) = -HA(1)
      HA(5) = -HA(2)
      HA(6) = HA(1)
      HA(7) = -HA(2)
      HA(8) = -HA(3)
      HA(9) = HA(3)
      HA(10) = HA(3)
      RETURN

  120 A = SQRT((X(3)-X(5))**2+ (X(4)-X(6))**2)
      B = (X(3)-X(5))/A
      C = (X(4)-X(6))/A
      HA(6) = -1.0D0/A + B*B/A
      HA(9) = B*C/A
      HA(10) = -1.0D0/A + C*C/A
      HA(13) = -HA(6)
      HA(14) = -HA(9)
      HA(15) = HA(6)
      HA(18) = -HA(9)
      HA(19) = -HA(10)
      HA(20) = HA(9)
      HA(21) = HA(10)
      RETURN

  130 A = SQRT((X(5)-X(1))**2+ (X(6)-X(2))**2)
      B = (X(5)-X(1))/A
      C = (X(6)-X(2))/A
      HA(1) = -1.0D0/A + B*B/A
      HA(2) = B*C/A
      HA(3) = -1.0D0/A + C*C/A
      HA(11) = -HA(1)
      HA(12) = -HA(2)
      HA(15) = HA(1)
      HA(16) = -HA(2)
      HA(17) = -HA(3)
      HA(20) = HA(2)
      HA(21) = HA(3)
      RETURN

  140 DO 150 I = 1,N* (N+1)/2
          HA(I) = 0.0D0
  150 CONTINUE
      DO 160 I = 1,N
          J = I* (I+1)/2
          HA(J) = -2.0D0*Y(KA)*Y(KA)*COS(Y(KA)*X(I))/1.5D1
  160 CONTINUE
      RETURN

  170 DO 180 I = 1,N* (N+1)/2
          HA(I) = 0.0D0
  180 CONTINUE
      K = 0
      DO 220 I = 1,5
          A = 0.0D0
          P = 0.0D0
          DO 190 J = 1,N
              A = A + Y(K+J)*X(J)** (1.0D0-Y(I+80))
              P = P + Y(K+J+40)*X(J)
  190     CONTINUE
          B = Y(K+KA)/ (X(KA)**Y(I+80)*A)
          C = B*Y(I+80)/X(KA)
          KK = 0
          DO 210 J = 1,N
              Q = Y(K+J)* (1.0D0-Y(I+80))/ (X(J)**Y(I+80)*A)
              DO 200 L = 1,J
                  R = Y(K+L)* (1.0D0-Y(I+80))/ (X(L)**Y(I+80)*A)
                  KK = KK + 1
                  HA(KK) = HA(KK) + B* (2.0D0*P*Q*R-Q*Y(K+L+40)-
     +                     R*Y(K+J+40))
                  IF (J.EQ.L) HA(KK) = HA(KK) + C*Q*P
                  IF (L.EQ.KA) HA(KK) = HA(KK) - C* (Y(K+J+40)-Q*P)
                  IF (J.EQ.KA) HA(KK) = HA(KK) - C* (Y(K+L+40)-R*P)
  200         CONTINUE
  210     CONTINUE
          KK = KA* (KA+1)/2
          Q = Y(K+KA)* (1.0D0-Y(I+80))/ (X(KA)**Y(I+80)*A)
          HA(KK) = HA(KK) + C*P* (1.0D0+Y(I+80))/X(KA)
          K = K + N
  220 CONTINUE
      RETURN

  230 DO 240 I = 1,N* (N+1)/2
          HA(I) = 0.0D0
  240 CONTINUE
      HA(1) = 2.0D0
      HA(2) = 1.0D0
      HA(3) = 2.0D0
      HA(6) = 2.0D0
      HA(10) = 8.0D0
      HA(15) = 2.0D0
      HA(21) = 4.0D0
      HA(28) = 1.0D1
      HA(36) = 1.4D1
      HA(45) = 4.0D0
      HA(55) = 2.0D0
      GO TO (300,250,260,270,280,290) KA

  250 HA(1) = HA(1) + 6.0D1
      HA(3) = HA(3) + 8.0D1
      HA(6) = HA(6) + 4.0D1
      RETURN

  260 HA(1) = HA(1) + 1.0D2
      HA(6) = HA(6) + 2.0D1
      RETURN

  270 HA(1) = HA(1) + 1.0D1
      HA(3) = HA(3) + 4.0D1
      HA(15) = HA(15) + 6.0D1
      RETURN

  280 HA(1) = HA(1) + 1.0D1
      HA(2) = HA(2) - 2.0D1
      HA(3) = HA(3) + 4.0D1
      RETURN

  290 HA(45) = HA(45) + 2.4D2
  300 RETURN

  310 DO 320 I = 1,N* (N+1)/2
          HA(I) = 0.0D0
  320 CONTINUE
      HA(1) = 2.0D0
      HA(2) = 1.0D0
      HA(3) = 2.0D0
      HA(6) = 2.0D0
      HA(10) = 8.0D0
      HA(15) = 2.0D0
      HA(21) = 4.0D0
      HA(28) = 1.0D1
      HA(36) = 1.4D1
      HA(45) = 4.0D0
      HA(55) = 2.0D0
      HA(66) = 2.0D0
      HA(78) = 2.0D1
      HA(91) = 1.0D1
      HA(105) = 8.0D0
      HA(120) = 5.4D1
      HA(136) = 1.2D1*X(16)**2
      HA(153) = 2.0D0
      HA(171) = 2.6D1
      HA(190) = 2.0D0
      HA(210) = 2.0D0
      GO TO (300,250,260,270,280,290,330,340,350,360,370,
     +       380,390,400) KA

  330 HA(1) = HA(1) + 2.0D1
      RETURN

  340 HA(91) = HA(91) + 1.0D2
      RETURN

  350 HA(91) = HA(91) + 6.0D1
      RETURN

  360 HA(1) = HA(1) + 2.8D2
      RETURN

  370 HA(3) = HA(3) + 3.0D2
      RETURN

  380 HA(1) = HA(1) + 1.0D2
      HA(153) = HA(153) + 10.8D2*X(17)**2
      RETURN

  390 HA(1) = HA(1) + 2.0D1
      RETURN

  400 HA(1) = HA(1) + 1.4D2
      HA(3) = HA(3) + 1.0D2
      HA(190) = HA(190) + 2.0D1
      RETURN

  410 DO 420 I = 1,N* (N+1)/2
          HA(I) = 0.0D0
  420 CONTINUE
      IF (MOD(KA,2).EQ.0) THEN
          I = (KA+2)/2
          J = I* (I+1)/2
          HA(J) = HA(J) + 2.0D0

      ELSE
          I = (KA+1)/2
          J = I* (I+1)/2
          HA(J) = HA(J) + 4.0D0
      END IF

      RETURN

  430 DO 440 I = 1,N* (N+1)/2
          HA(I) = 0.0D0
  440 CONTINUE
      K = (N-1)/2
      KK = K* (K+1)/2
      A = Y(KA)**2
      L = 0
      LL = KK
      S = 1.0D0
      IF (KA.GT.62 .AND. KA.LT.124) S = -S
      DO 450 I = 1,K
          L = L + I
          B = X(I+K)**2 - A
          C = A*X(I)**2
          P = B*B + C
          Q = B*B - C
          R = S*P*P
          HA(L) = A*Q/R
          HA(L+LL) = -4.0D0*A*B*X(I)*X(I+K)/R
          LL = LL + K
          HA(L+LL) = 2.0D0*S*B/P - 4.0D0*Q*X(I+K)**2/R
  450 CONTINUE
      RETURN

  460 P = 5.0D2
      A = 0.99D0
      DO 470 I = 1,N* (N+1)/2
          HA(I) = 0.0D0
  470 CONTINUE
      HA(25) = -0.63D-1
      GO TO (480,490,500,510,520,530,540,550,560) KA

  480 CONTINUE
      RETURN

  490 HA(29) = P* (1.3167D-1-1.334D-2*X(8))
      HA(36) = -P*X(1)*1.334D-2
      RETURN

  500 HA(29) = -P* (1.3167D-1-1.334D-2*X(8))
      HA(36) = P*X(1)*1.334D-2
      RETURN

  510 HA(36) = -0.76D-1*P
      RETURN

  520 HA(36) = 0.76D-1*P
      RETURN

  530 C = (X(4)*X(9)+1.0D3*X(3))**3
      Q = (X(4)*X(9)-1.0D3*X(3))
      HA(6) = -1.96D8*P*X(4)*X(9)/C
      HA(9) = -9.80D4*P*X(9)*Q/C
      HA(10) = 1.96D5*P*X(3)*X(9)**2/C
      HA(39) = -9.80D4*P*X(4)*Q/C
      HA(40) = 9.80D4*P*X(3)*Q/C
      HA(45) = 1.96D5*P*X(3)*X(4)**2/C
      RETURN

  540 C = (X(4)*X(9)+1.0D3*X(3))**3
      Q = (X(4)*X(9)-1.0D3*X(3))
      HA(6) = 1.96D8*P*X(4)*X(9)/C
      HA(9) = 9.80D4*P*X(9)*Q/C
      HA(10) = -1.96D5*P*X(3)*X(9)**2/C
      HA(39) = 9.80D4*P*X(4)*Q/C
      HA(40) = -9.80D4*P*X(3)*Q/C
      HA(45) = -1.96D5*P*X(3)*X(4)**2/C
      RETURN

  550 HA(1) = 2.0D0*P* (X(2)+X(5))/X(1)**3
      HA(2) = -P/X(1)**2
      HA(11) = -P/X(1)**2
      RETURN

  560 HA(1) = -2.0D0*P* (X(2)+X(5))/X(1)**3
      HA(2) = P/X(1)**2
      HA(11) = P/X(1)**2
      RETURN

  570 P = 1.0D5
      DO 580 I = 1,N* (N+1)/2
          HA(I) = 0.0D0
  580 CONTINUE
      HA(13) = Y(6)
      HA(16) = Y(2)
      GO TO (730,590,600,610,620,630,640,650,660,
     +       670,680,690,700) KA

  590 HA(1) = 2.0D0*P*Y(8)*X(3)/X(1)**3
      HA(4) = -P*Y(8)/X(1)**2
      HA(21) = 2.0D0*P*Y(7)
      RETURN

  600 HA(4) = -P* (Y(10)+Y(11)*X(6)+Y(12)*X(6)**2)/X(3)**2
      HA(6) = 2.0D0*P* (Y(10)+Y(11)*X(6)+Y(12)*X(6)**2)*X(1)/X(3)**3
      HA(16) = HA(16) + P* (Y(11)+2.0D0*Y(12)*X(6))/X(3)
      HA(18) = -P* (Y(11)+2.0D0*Y(12)*X(6))*X(1)/X(3)**2
      HA(21) = 2.0D0*P*Y(12)*X(1)/X(3)
      RETURN

  610 HA(21) = 2.0D0*P*Y(13)
      RETURN

  620 HA(14) = -P*Y(19)/X(5)**2
      HA(15) = 2.0D0*P* (Y(17)+Y(18)*X(6)+Y(19)*X(4)+Y(20)*X(6)**2)/
     +         X(5)**3
      HA(20) = -P* (Y(18)+2.0D0*Y(20)*X(6))/X(5)**2
      HA(21) = 2.0D0*P*Y(20)/X(5)
      RETURN

  630 HA(5) = -P* (Y(22)/X(4)+Y(23))/X(3)**2
      HA(6) = 2.0D0*P* (Y(22)/X(4)+Y(23))*X(2)/X(3)**3
      HA(8) = -P*Y(22)/ (X(3)*X(4)**2)
      HA(9) = P*Y(22)*X(2)/ (X(3)*X(4))**2
      HA(10) = 2.0D0*P*Y(22)*X(2)/ (X(3)*X(4)**3)
      RETURN

  640 HA(5) = -P* (Y(25)+Y(26)/X(4))/ (X(3)**2*X(7))
      HA(6) = 2.0D0*P* (Y(25)+Y(26)/X(4))*X(2)/ (X(3)**3*X(7))
      HA(8) = -P*Y(26)/ (X(3)*X(4)**2*X(7))
      HA(9) = P*Y(26)*X(2)/ ((X(3)*X(4))**2*X(7))
      HA(10) = 2.0D0*P*Y(26)*X(2)/ (X(3)*X(4)**3*X(7))
      HA(23) = -P* (Y(25)+Y(26)/X(4))/ (X(3)*X(7)**2)
      HA(24) = P* (Y(25)+Y(26)/X(4))*X(2)/ (X(3)*X(7))**2
      HA(25) = P*Y(26)*X(2)/ (X(3)* (X(4)*X(7))**2)
      HA(28) = 2.0D0*P* (Y(24)+ (Y(25)+Y(26)/X(4))*X(2)/X(3))/X(7)**3
      RETURN

  650 HA(15) = 2.0D0*P* (Y(27)+Y(28)*X(7))/X(5)**3
      HA(26) = -P*Y(28)/X(5)**2
      RETURN

  660 HA(4) = -P*Y(33)/X(3)**2
      HA(6) = 2.0D0*P* (Y(33)*X(1)+Y(34))/X(3)**3
      RETURN

  670 HA(5) = -P* (Y(35)/X(4)+Y(36))/X(3)**2
      HA(6) = 2.0D0*P* (Y(35)/X(4)+Y(36))*X(2)/X(3)**3
      HA(8) = -P*Y(35)/ (X(3)*X(4)**2)
      HA(9) = P*Y(35)*X(2)/ (X(3)*X(4))**2
      HA(10) = 2.0D0*P*Y(35)*X(2)/ (X(3)*X(4)**3)
      RETURN

  680 HA(3) = 2.0D0*P*Y(38)*X(3)*X(4)/X(2)**3
      HA(5) = -P*Y(38)*X(4)/X(2)**2
      HA(8) = -P*Y(38)*X(3)/X(2)**2
      HA(9) = P*Y(38)/X(2)
      RETURN

  690 HA(16) = HA(16) + P*Y(39)
      RETURN

  700 HA(1) = 2.0D0*P* (Y(42)*X(3)+Y(43))/X(1)**3
      HA(4) = -P*Y(42)/X(1)**2
      RETURN

  710 P = 1.0D5
      DO 720 I = 1,N* (N+1)/2
          HA(I) = 0.0D0
  720 CONTINUE
      GO TO (730,740,750,760) KA

  730 CONTINUE
      RETURN

  740 HA(1) = 2.0D0*P* (833.33252D0*X(4)-83333.333D0)/ (X(1)**3*X(6))
      HA(7) = -P*833.33252D0/ (X(1)**2*X(6))
      HA(16) = P* (833.33252D0*X(4)-83333.333D0)/ (X(1)**2*X(6)**2)
      HA(19) = -P*833.33252D0/ (X(1)*X(6)**2)
      HA(21) = 2.0D0*P* (833.33252D0*X(4)/X(1)+1.0D2-83333.333D0/X(1))/
     +         X(6)**3
      RETURN

  750 HA(3) = P*2.50D3* (X(5)-X(4))/ (X(2)**3*X(7))
      HA(8) = P*1.25D3/ (X(2)**2*X(7))
      HA(12) = -P*1.25D3/ (X(2)**2*X(7))
      HA(23) = P*1.25D3* (X(5)-X(4))/ (X(2)**2*X(7)**2)
      HA(25) = -P* (1.0D0-1.25D3/X(2))/X(7)**2
      HA(26) = -P*1.25D3/ (X(2)*X(7)**2)
      HA(28) = 2.0D0*P* (1.25D3* (X(5)-X(4))/X(2)+X(4))/X(7)**3
      RETURN

  760 HA(6) = 2.0D0*P* (1.25D6-2.5D3*X(5))/ (X(3)**3*X(8))
      HA(13) = P*2.5D3/ (X(3)**2*X(8))
      HA(31) = P* (1.25D6-2.5D3*X(5))/ (X(3)**2*X(8)**2)
      HA(33) = -P* (1.00D0-2.5D3/X(3))/X(8)**2
      HA(36) = 2.0D0*P* ((1.25D6-2.5D3*X(5))/X(3)+X(5))/X(8)**3
      RETURN

  770 P = 2.0D3
      DO 780 I = 1,N* (N+1)/2
          HA(I) = 0.0D0
  780 CONTINUE
      HA(67) = -1.231060D0
      HA(80) = -1.231060D0
      HA(94) = -1.231060D0
      HA(109) = -1.231060D0
      HA(125) = -1.231060D0
      GO TO (730,790,800,810,820,830,840,850,860,870,
     +       880,890,900,910,920,930,940,950,960) KA

  790 HA(1) = -P*1.95D-2/X(6)
      HA(16) = -P* (3.475D-2-1.95D-2*X(1))/X(6)**2
      HA(21) = 2.0D0*P*X(1)* (3.475D-2-9.75D-3*X(1))/X(6)**3
      RETURN

  800 HA(3) = -P*1.95D-2/X(7)
      HA(23) = -P* (3.475D-2-1.95D-2*X(2))/X(7)**2
      HA(28) = 2.0D0*P*X(2)* (3.475D-2-9.75D-3*X(2))/X(7)**3
      RETURN

  810 HA(6) = -P*1.95D-2/X(8)
      HA(31) = -P* (3.475D-2-1.95D-2*X(3))/X(8)**2
      HA(36) = 2.0D0*P*X(3)* (3.475D-2-9.75D-3*X(3))/X(8)**3
      RETURN

  820 HA(10) = -P*1.95D-2/X(9)
      HA(40) = -P* (3.475D-2-1.95D-2*X(4))/X(9)**2
      HA(45) = 2.0D0*P*X(4)* (3.475D-2-9.75D-3*X(4))/X(9)**3
      RETURN

  830 HA(15) = -P*1.95D-2/X(10)
      HA(50) = -P* (3.475D-2-1.95D-2*X(5))/X(10)**2
      HA(55) = 2.0D0*P*X(5)* (3.475D-2-9.75D-3*X(5))/X(10)**3
      RETURN

  840 HA(22) = -P*X(12)/ (X(11)*X(7)**2)
      HA(27) = -P* (1.0D0-X(12)/X(11))/X(7)**2
      HA(28) = 2.0D0*P* (X(6)+ (X(1)-X(6))*X(12)/X(11))/X(7)**3
      HA(56) = -P*X(12)/ (X(7)*X(11)**2)
      HA(61) = P*X(12)/ (X(7)*X(11)**2)
      HA(62) = P* (X(1)-X(6))*X(12)/ (X(11)*X(7))**2
      HA(66) = 2.0D0*P* (X(1)-X(6))*X(12)/ (X(7)*X(11)**3)
      HA(67) = HA(67) + P/ (X(11)*X(7))
      HA(72) = -P/ (X(11)*X(7))
      HA(73) = -P* (X(1)-X(6))/ (X(11)*X(7)**2)
      HA(77) = -P* (X(1)-X(6))/ (X(7)*X(11)**2)
      RETURN

  850 HA(29) = 2.0D-3*P*X(12)/X(8)**2
      HA(30) = -2.0D-3*P*X(13)/X(8)**2
      HA(35) = -P* (1.0D0+2.0D-3*X(12))/X(8)**2
      HA(36) = 2.0D0*P* (X(7)+2.0D-3* ((X(7)-X(1))*X(12)+X(2)*X(13)))/
     +         X(8)**3
      HA(67) = HA(67) - 2.0D-3*P/X(8)
      HA(73) = 2.0D-3*P/X(8)
      HA(74) = -2.0D-3*P* (X(7)-X(1))/X(8)**2
      HA(80) = HA(80) + 2.0D-3*P/X(8)
      HA(86) = -2.0D-3*P*X(2)/X(8)**2
      RETURN

  860 HA(80) = HA(80) - P*2.0D-3
      HA(86) = P*2.0D-3
      HA(94) = HA(94) + P*2.0D-3
      HA(100) = -P*2.0D-3
      RETURN

  870 HA(6) = 2.0D0*P* (X(9)+ ((X(4)-X(8))*X(15)+5.0D2* (X(10)-X(9)))/
     +        X(14))/X(3)**3
      HA(9) = -P*X(15)/ (X(14)*X(3)**2)
      HA(31) = P*X(15)/ (X(14)*X(3)**2)
      HA(39) = -P* (1.0D0-5.0D2/X(14))/X(3)**2
      HA(48) = -5.0D2*P/ (X(14)*X(3)**2)
      HA(94) = HA(94) + P* ((X(4)-X(8))*X(15)+5.0D2* (X(10)-X(9)))/
     +         (X(14)*X(3))**2
      HA(95) = -P*X(15)/ (X(3)*X(14)**2)
      HA(99) = P*X(15)/ (X(3)*X(14)**2)
      HA(100) = 5.0D2*P/ (X(3)*X(14)**2)
      HA(101) = -5.0D2*P/ (X(3)*X(14)**2)
      HA(105) = 2.0D0*P* ((X(4)-X(8))*X(15)+5.0D2* (X(10)-X(9)))/
     +          (X(3)*X(14)**3)
      HA(108) = -P* (X(4)-X(8))/ (X(14)*X(3)**2)
      HA(109) = HA(109) + P/ (X(3)*X(14))
      HA(113) = -P/ (X(3)*X(14))
      HA(119) = -P* (X(4)-X(8))/ (X(3)*X(14)**2)
      RETURN

  880 HA(10) = 2.0D0*P* (X(10)+ (X(5)*X(16)-5.0D2*X(10))/X(15))/X(4)**3
      HA(14) = -P*X(16)/ (X(15)*X(4)**2)
      HA(49) = -P* (1.0D0-5.0D2/X(15))/X(4)**2
      HA(109) = HA(109) + P* (X(5)*X(16)-5.0D2*X(10))/ (X(15)*X(4))**2
      HA(110) = -P*X(16)/ (X(4)*X(15)**2)
      HA(115) = 5.0D2*P/ (X(4)*X(15)**2)
      HA(120) = 2.0D0*P* (5.0D2-X(16)+ (X(5)*X(16)-5.0D2*X(10))/X(4))/
     +          X(15)**3
      HA(124) = -P*X(5)/ (X(15)*X(4)**2)
      HA(125) = HA(125) + P/ (X(4)*X(15))
      HA(135) = -P* (X(5)/X(4)-1.0D0)/X(15)**2
      RETURN

  890 HA(10) = 2.0D0*P* (9.0D-1-2.0D-3*X(5)*X(16))/X(4)**3
      HA(14) = P* (2.0D-3*X(16))/X(4)**2
      HA(124) = P* (2.0D-3*X(5))/X(4)**2
      HA(125) = HA(125) - P*2.0D-3/X(4)
      RETURN

  900 HA(66) = 2.0D0*P*X(12)/X(11)**3
      HA(77) = -P/X(11)**2
      RETURN

  910 HA(14) = -P/X(5)**2
      HA(15) = 2.0D0*P*X(4)/X(5)**3
      RETURN

  920 HA(9) = -P/X(4)**2
      HA(10) = 2.0D0*P*X(3)/X(4)**3
      RETURN

  930 HA(5) = -P/X(3)**2
      HA(6) = 2.0D0*P*X(2)/X(3)**3
      RETURN

  940 HA(2) = -P/X(2)**2
      HA(3) = 2.0D0*P*X(1)/X(2)**3
      RETURN

  950 HA(54) = -P/X(10)**2
      HA(55) = 2.0D0*P*X(9)/X(10)**3
      RETURN

  960 HA(44) = -P/X(9)**2
      HA(45) = 2.0D0*P*X(8)/X(9)**3
      RETURN

      END
* SUBROUTINE TYTIM1                MS DOS                     91/12/01
C PORTABILITY : MS DOS / MS FORTRAN v.5.0
C 91/12/01 SI : ORIGINAL VERSION
*
* PURPOSE :
*  GET TIME IN 100TH OF SEC.
*
      SUBROUTINE TYTIM1(ITIME)
C     .. Scalar Arguments ..
      INTEGER ITIME
C     ..
C     .. Local Scalars ..
C     INTEGER*2 I100TH,IHR,IMIN,ISEC
C     ..
C     .. External Subroutines ..
C     EXTERNAL GETTIM
C     ..
C     CALL GETTIM(IHR,IMIN,ISEC,I100TH)
C     ITIME = 100* (IHR*60*60+IMIN*60+ISEC) + I100TH
      ITIME = 0
      END
* SUBROUTINE TYTIM2                ALL SYSTEMS                91/12/01
C PORTABILITY : ALL SYSTEMS
C 91/12/01 SI : ORIGINAL VERSION
*
* PURPOSE :
*  PRINT TIME ELAPSED.
*
      SUBROUTINE TYTIM2(ITIME)
C     .. Scalar Arguments ..
      INTEGER ITIME
C     ..
C     .. Local Scalars ..
      INTEGER IHR,IMIN,ISEC,IT
C     ..
C     .. External Subroutines ..
      EXTERNAL TYTIM1
C     ..
C     CALL TYTIM1(IT)
C     IT = IT - ITIME
C     IHR = IT/ (60*60*100)
C     IT = IT - IHR*60*60*100
C     IMIN = IT/ (60*100)
C     IT = IT - IMIN*60*100
C     ISEC = IT/100
C     IT = IT - ISEC*100
C     WRITE (6,FMT=9000) IHR,IMIN,ISEC,IT

C9000 FORMAT (' TIME=',I2,':',I2.2,':',I2.2,'.',I2.2)
      END
