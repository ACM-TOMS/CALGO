      SUBROUTINE FILL0(BITMAP, N)                                       FIL   10
C
C     FILL THE FIRST N BITS OF BITMAP WITH ZEROES.
C
      INTEGER BITMAP(1), N
C
      DATA NBPW /35/
C     NBPW IS THE MINIMUM NUMBER OF SIGNIFICANT BITS PER WORD USED
C     BY INTEGER ARITHMETIC.  THIS IS USUALLY ONE LESS THAN THE
C     ACTUAL NUMBER OF BITS PER WORD, BUT AN IMPORTANT EXCEPTION IS
C     THE CDC-6000 SERIES OF MACHINES, WHERE NBPW SHOULD BE 48.
C
      LOOP = N/NBPW
      NBLW = MOD(N,NBPW)
      IF (LOOP.EQ.0) GO TO 20
      DO 10 I=1,LOOP
        BITMAP(I) = 0
   10 CONTINUE
   20 IF (NBLW.NE.0) BITMAP(LOOP+1) = MOD(BITMAP(LOOP+1),2**(NBPW-NBLW)
     *  )
      RETURN
      END
      SUBROUTINE MARK1(BITMAP, N)                                       MAR   10
C
C     PUT A ONE IN THE NTH BIT OF BITMAP.
C
      INTEGER BITMAP(1), N
C
      DATA NBPW /35/
C     NBPW IS THE MINIMUM NUMBER OF SIGNIFICANT BITS PER WORD USED
C     BY INTEGER ARITHMETIC.  THIS IS USUALLY ONE LESS THAN THE
C     ACTUAL NUMBER OF BITS PER WORD, BUT AN IMPORTANT EXCEPTION IS
C     THE CDC-6000 SERIES OF MACHINES, WHERE NBPW SHOULD BE 48.
C
      NWORD = (N-1)/NBPW
      NBIT = MOD(N-1,NBPW)
      I = 2**(NBPW-NBIT-1)
      BITMAP(NWORD+1) = BITMAP(NWORD+1) + I*(1-MOD(BITMAP(NWORD+1)/I,2)
     *  )
      RETURN
      END
      FUNCTION IGET(BITMAP, N)                                          IGE   10
C
C     IGET=0 IF THE NTH BIT OF BITMAP IS ZERO, ELSE IGET IS ONE.
C
      INTEGER BITMAP(1), N
C
      DATA NBPW /35/
C     NBPW IS THE MINIMUM NUMBER OF SIGNIFICANT BITS PER WORD USED
C     BY INTEGER ARITHMETIC.  THIS IS USUALLY ONE LESS THAN THE
C     ACTUAL NUMBER OF BITS PER WORD, BUT AN IMPORTANT EXCEPTION IS
C     THE CDC-6000 SERIES OF MACHINES, WHERE NBPW SHOULD BE 48.
C
      NWORD = (N-1)/NBPW
      NBIT = MOD(N-1,NBPW)
      IGET = MOD(BITMAP(NWORD+1)/2**(NBPW-NBIT-1),2)
      RETURN
      END
      SUBROUTINE GCONTR(Z, NRZ, NX, NY, CV, NCV, ZMAX, BITMAP, DRAW)    GCO   10
C
C     THIS SUBROUTINE DRAWS A CONTOUR THROUGH EQUAL VALUES OF AN ARRAY.
C
C     *****     FORMAL ARGUMENTS     ***********************************
C
C     Z IS THE ARRAY FOR WHICH CONTOURS ARE TO BE DRAWN.  THE ELEMENTS
C     OF Z ARE ASSUMED TO LIE UPON THE NODES OF A TOPOLOGICALLY
C     RECTANGULAR COORDINATE SYSTEM - E.G. CARTESIAN, POLAR (EXCEPT
C     THE ORIGIN), ETC.
C
C     NRZ IS THE NUMBER OF ROWS DECLARED FOR Z IN THE CALLING PROGRAM.
C
C     NX IS THE LIMIT FOR THE FIRST SUBSCRIPT OF Z.
C
C     NY IS THE LIMIT FOR THE SECOND SUBSCRIPT OF Z.
C
C     CV ARE THE VALUES OF THE CONTOURS TO BE DRAWN.
C
C     NCV IS THE NUMBER OF CONTOUR VALUES IN CV.
C
C     ZMAX IS THE MAXIMUM VALUE OF Z FOR CONSIDERATION.  A VALUE OF
C     Z(I,J) GREATER THAN ZMAX IS A SIGNAL THAT THAT POINT AND THE
C     GRID LINE SEGMENTS RADIATING FROM THAT POINT TO IT'S NEIGHBORS
C     ARE TO BE EXCLUDED FROM CONTOURING.
C
C     BITMAP IS A WORK AREA LARGE ENOUGH TO HOLD 2*NX*NY*NCV BITS.  IT
C     IS ACCESSED BY LOW-LEVEL ROUTINES, WHICH ARE DESCRIBED BELOW.
C     LET J BE THE NUMBER OF USEFUL BITS IN EACH WORD OF BITMAP,
C     AS DETERMINED BY THE USER MACHINE AND IMPLEMENTATION OF
C     THE BITMAP MANIPULATION SUBPROGRAMS DESCRIBED BELOW.  THEN
C     THE NUMBER OF WORDS REQUIRED FOR THE BITMAP IS THE FLOOR OF
C         (2*NX*NY*NCV+J-1)/J.
C
C     DRAW IS A USER-PROVIDED SUBROUTINE USED TO DRAW CONTOURS.
C     THE CALLING SEQUENCE FOR DRAW IS:
C
C         CALL DRAW (X,Y,IFLAG)
C         LET NX = INTEGER PART OF X, FX = FRACTIONAL PART OF X.
C         THEN X SHOULD BE INTERPRETED SUCH THAT INCREASES IN NX
C         CORRESPOND TO INCREASES IN THE FIRST SUBSCRIPT OF Z, AND
C         FX IS THE FRACTIONAL DISTANCE FROM THE ABSCISSA CORRESPONDING
C         TO NX TO THE ABSCISSA CORRESPONDING TO NX+1,
C         AND Y SHOULD BE INTERPRETED SIMILARLY FOR THE SECOND
C         SUBSCRIPT OF Z.
C         THE LOW-ORDER DIGIT OF IFLAG WILL HAVE ONE OF THE VALUES:
C             1 - CONTINUE A CONTOUR,
C             2 - START A CONTOUR AT A BOUNDARY,
C             3 - START A CONTOUR NOT AT A BOUNDARY,
C             4 - FINISH A CONTOUR AT A BOUNDARY,
C             5 - FINISH A CLOSED CONTOUR (NOT AT A BOUNDARY).
C                 NOTE THAT REQUESTS 1, 4 AND 5 ARE FOR PEN-DOWN
C                 MOVES, AND THAT REQUESTS 2 AND 3 ARE FOR PEN-UP
C                 MOVES.
C             6 - SET X AND Y TO THE APPROXIMATE 'PEN' POSITION, USING
C                 THE NOTATION DISCUSSED ABOVE.  THIS CALL MAY BE
C                 IGNORED, THE RESULT BEING THAT THE 'PEN' POSITION
C                 IS TAKEN TO CORRESPOND TO Z(1,1).
C         IFLAG/10 IS THE CONTOUR NUMBER.
C
C     *****     EXTERNAL SUBPROGRAMS     *******************************
C
C     DRAW IS THE USER-SUPPLIED LINE DRAWING SUBPROGRAM DESCRIBED ABOVE.
C     DRAW MAY BE SENSITIVE TO THE HOST COMPUTER AND TO THE PLOT DEVICE.
C     FILL0 IS USED TO FILL A BITMAP WITH ZEROES.  CALL FILL0 (BITMAP,N)
C     FILLS THE FIRST N BITS OF BITMAP WITH ZEROES.
C     MARK1 IS USED TO PLACE A 1 IN A SPECIFIC BIT OF THE BITMAP.
C     CALL MARK1 (BITMAP,N) PUTS A 1 IN THE NTH BIT OF THE BITMAP.
C     IGET IS USED TO DETERMINE THE SETTING OF A PARTICULAR BIT IN THE
C     BITMAP.  I=IGET(BITMAP,N) SETS I TO ZERO IF THE NTH BIT OF THE
C     BITMAP IS ZERO, AND SETS I TO ONE IF THE NTH BIT IS ONE.
C     FILL0, MARK1 AND IGET ARE MACHINE SENSITIVE.
C
C     ******************************************************************
C
      REAL Z(NRZ,1), CV(1)
      INTEGER BITMAP(1)
      INTEGER L1(4), L2(4), IJ(2)
C
C     L1 AND L2 CONTAIN LIMITS USED DURING THE SPIRAL SEARCH FOR THE
C     BEGINNING OF A CONTOUR.
C     IJ STORES SUBCRIPTS USED DURING THE SPIRAL SEARCH.
C
      INTEGER I1(2), I2(2), I3(6)
C
C     I1, I2 AND I3 ARE USED FOR SUBSCRIPT COMPUTATIONS DURING THE
C     EXAMINATION OF LINES FROM Z(I,J) TO IT'S NEIGHBORS.
C
      REAL XINT(4)
C
C     XINT IS USED TO MARK INTERSECTIONS OF THE CONTOUR UNDER
C     CONSIDERATION WITH THE EDGES OF THE CELL BEING EXAMINED.
C
      REAL XY(2)
C
C     XY IS USED TO COMPUTE COORDINATES FOR THE DRAW SUBROUTINE.
C
      EQUIVALENCE (L2(1),IMAX), (L2(2),JMAX), (L2(3),IMIN),
     *  (L2(4),JMIN)
      EQUIVALENCE (IJ(1),I), (IJ(2),J)
      EQUIVALENCE (XY(1),X), (XY(2),Y)
C
      DATA L1(3) /-1/, L1(4) /-1/
      DATA I1 /1,0/, I2 /1,-1/, I3 /1,0,0,1,1,0/
C
      L1(1) = NX
      L1(2) = NY
      DMAX = ZMAX
C
C     SET THE CURRENT PEN POSITION.  THE DEFAULT POSITION CORRESPONDS
C     TO Z(1,1).
C
      X = 1.0
      Y = 1.0
      CALL DRAW(X, Y, 6)
      ICUR = MAX0(1,MIN0(INT(X),NX))
      JCUR = MAX0(1,MIN0(INT(Y),NY))
C
C     CLEAR THE BITMAP
C
      CALL FILL0(BITMAP, 2*NX*NY*NCV)
C
C     SEARCH ALONG A RECTANGULAR SPIRAL PATH FOR A LINE SEGMENT HAVING
C     THE FOLLOWING PROPERTIES:
C          1.  THE END POINTS ARE NOT EXCLUDED,
C          2.  NO MARK HAS BEEN RECORDED FOR THE SEGMENT,
C          3.  THE VALUES OF Z AT THE ENDS OF THE SEGMENT ARE SUCH THAT
C              ONE Z IS LESS THAN THE CURRENT CONTOUR VALUE, AND THE
C              OTHER IS GREATER THAN OR EQUAL TO THE CURRENT CONTOUR
C              VALUE.
C
C     SEARCH ALL BOUNDARIES FIRST, THEN SEARCH INTERIOR LINE SEGMENTS.
C     NOTE THAT THE INTERIOR LINE SEGMENTS NEAR EXCLUDED POINTS MAY BE
C     BOUNDARIES.
C
      IBKEY = 0
   10 I = ICUR
      J = JCUR
   20 IMAX = I
      IMIN = -I
      JMAX = J
      JMIN = -J
      IDIR = 0
C     DIRECTION ZERO IS +I, 1 IS +J, 2 IS -I, 3 IS -J.
   30 NXIDIR = IDIR + 1
      K = NXIDIR
      IF (NXIDIR.GT.3) NXIDIR = 0
   40 I = IABS(I)
      J = IABS(J)
      IF (Z(I,J).GT.DMAX) GO TO 140
      L = 1
C     L=1 MEANS HORIZONTAL LINE, L=2 MEANS VERTICAL LINE.
   50 IF (IJ(L).GE.L1(L)) GO TO 130
      II = I + I1(L)
      JJ = J + I1(3-L)
      IF (Z(II,JJ).GT.DMAX) GO TO 130
      ASSIGN 100 TO JUMP
C     THE NEXT 15 STATEMENTS (OR SO) DETECT BOUNDARIES.
   60 IX = 1
      IF (IJ(3-L).EQ.1) GO TO 80
      II = I - I1(3-L)
      JJ = J - I1(L)
      IF (Z(II,JJ).GT.DMAX) GO TO 70
      II = I + I2(L)
      JJ = J + I2(3-L)
      IF (Z(II,JJ).LT.DMAX) IX = 0
   70 IF (IJ(3-L).GE.L1(3-L)) GO TO 90
   80 II = I + I1(3-L)
      JJ = J + I1(L)
      IF (Z(II,JJ).GT.DMAX) GO TO 90
      IF (Z(I+1,J+1).LT.DMAX) GO TO JUMP, (100, 280)
   90 IX = IX + 2
      GO TO JUMP, (100, 280)
  100 IF (IX.EQ.3) GO TO 130
      IF (IX+IBKEY.EQ.0) GO TO 130
C     NOW DETERMINE WHETHER THE LINE SEGMENT IS CROSSED BY THE CONTOUR.
      II = I + I1(L)
      JJ = J + I1(3-L)
      Z1 = Z(I,J)
      Z2 = Z(II,JJ)
      DO 120 ICV=1,NCV
        IF (IGET(BITMAP,2*(NX*(NY*(ICV-1)+J-1)+I-1)+L).NE.0) GO TO 120
        IF (CV(ICV).LE.AMIN1(Z1,Z2)) GO TO 110
        IF (CV(ICV).LE.AMAX1(Z1,Z2)) GO TO 190
  110   CALL MARK1(BITMAP, 2*(NX*(NY*(ICV-1)+J-1)+I-1)+L)
  120 CONTINUE
  130 L = L + 1
      IF (L.LE.2) GO TO 50
  140 L = MOD(IDIR,2) + 1
      IJ(L) = ISIGN(IJ(L),L1(K))
C
C     LINES FROM Z(I,J) TO Z(I+1,J) AND Z(I,J+1) ARE NOT SATISFACTORY.
C     CONTINUE THE SPIRAL.
C
  150 IF (IJ(L).GE.L1(K)) GO TO 170
      IJ(L) = IJ(L) + 1
      IF (IJ(L).GT.L2(K)) GO TO 160
      GO TO 40
  160 L2(K) = IJ(L)
      IDIR = NXIDIR
      GO TO 30
  170 IF (IDIR.EQ.NXIDIR) GO TO 180
      NXIDIR = NXIDIR + 1
      IJ(L) = L1(K)
      K = NXIDIR
      L = 3 - L
      IJ(L) = L2(K)
      IF (NXIDIR.GT.3) NXIDIR = 0
      GO TO 150
  180 IF (IBKEY.NE.0) RETURN
      IBKEY = 1
      GO TO 10
C
C     AN ACCEPTABLE LINE SEGMENT HAS BEEN FOUND.
C     FOLLOW THE CONTOUR UNTIL IT EITHER HITS A BOUNDARY OR CLOSES.
C
  190 IEDGE = L
      CVAL = CV(ICV)
      IF (IX.NE.1) IEDGE = IEDGE + 2
      IFLAG = 2 + IBKEY
      XINT(IEDGE) = (CVAL-Z1)/(Z2-Z1)
  200 XY(L) = FLOAT(IJ(L)) + XINT(IEDGE)
      XY(3-L) = FLOAT(IJ(3-L))
      CALL MARK1(BITMAP, 2*(NX*(NY*(ICV-1)+J-1)+I-1)+L)
      CALL DRAW(X, Y, IFLAG+10*ICV)
      IF (IFLAG.LT.4) GO TO 210
      ICUR = I
      JCUR = J
      GO TO 20
C
C     CONTINUE A CONTOUR.  THE EDGES ARE NUMBERED CLOCKWISE WITH
C     THE BOTTOM EDGE BEING EDGE NUMBER ONE.
C
  210 NI = 1
      IF (IEDGE.LT.3) GO TO 220
      I = I - I3(IEDGE)
      J = J - I3(IEDGE+2)
  220 DO 250 K=1,4
        IF (K.EQ.IEDGE) GO TO 250
        II = I + I3(K)
        JJ = J + I3(K+1)
        Z1 = Z(II,JJ)
        II = I + I3(K+1)
        JJ = J + I3(K+2)
        Z2 = Z(II,JJ)
        IF (CVAL.LE.AMIN1(Z1,Z2)) GO TO 250
        IF (CVAL.GT.AMAX1(Z1,Z2)) GO TO 250
        IF (K.EQ.1) GO TO 230
        IF (K.NE.4) GO TO 240
  230   ZZ = Z1
        Z1 = Z2
        Z2 = ZZ
  240   XINT(K) = (CVAL-Z1)/(Z2-Z1)
        NI = NI + 1
        KS = K
  250 CONTINUE
      IF (NI.EQ.2) GO TO 260
C
C     THE CONTOUR CROSSES ALL FOUR EDGES OF THE CELL BEING EXAMINED.
C     CHOOSE THE LINES TOP-TO-LEFT AND BOTTOM-TO-RIGHT IF THE
C     INTERPOLATION POINT ON THE TOP EDGE IS LESS THAN THE INTERPOLATION
C     POINT ON THE BOTTOM EDGE.  OTHERWISE, CHOOSE THE OTHER PAIR.  THIS
C     METHOD PRODUCES THE SAME RESULTS IF THE AXES ARE REVERSED.  THE
C     CONTOUR MAY CLOSE AT ANY EDGE, BUT MUST NOT CROSS ITSELF INSIDE
C     ANY CELL.
C
      KS = 5 - IEDGE
      IF (XINT(3).LT.XINT(1)) GO TO 260
      KS = 3 - IEDGE
      IF (KS.LE.0) KS = KS + 4
C
C     DETERMINE WHETHER THE CONTOUR WILL CLOSE OR RUN INTO A BOUNDARY
C     AT EDGE KS OF THE CURRENT CELL.
C
  260 L = KS
      IFLAG = 1
      ASSIGN 280 TO JUMP
      IF (KS.LT.3) GO TO 270
      I = I + I3(KS)
      J = J + I3(KS+2)
      L = KS - 2
  270 IF (IGET(BITMAP,2*(NX*(NY*(ICV-1)+J-1)+I-1)+L).EQ.0) GO TO 60
      IFLAG = 5
      GO TO 290
  280 IF (IX.NE.0) IFLAG = 4
  290 IEDGE = KS + 2
      IF (IEDGE.GT.4) IEDGE = IEDGE - 4
      XINT(IEDGE) = XINT(KS)
      GO TO 200
C
      END
      DIMENSION Z(51,51), C(10), WORK(1680)                             MAN   10
C     DIMENSION OF WORK IS LARGE ENOUGH TO CONTAIN                      MAN   20
C     2*(DIMENSION OF C)*(TOTAL DIMENSION OF Z) USEFUL BITS.  SEE THE   MAN   30
C     BITMAP ROUTINES ACCESSED BY GCONTR.                               MAN   40
      REAL MU                                                           MAN   50
      EXTERNAL DRAW                                                     MAN   60
      COMMON /CUR/ XCUR, YCUR                                           MAN   70
      DATA C(1), C(2), C(3), C(4), C(5) /3.05,3.2,3.5,3.50135,3.6/      MAN   80
      DATA C(6), C(7), C(8), C(9), C(10) /3.766413,4.0,4.130149,5.0,    MAN   90
     *  10.0/                                                           MAN  100
      DATA NX /51/, NY /51/, NF /10/                                    MAN  110
      DATA XMIN /-2.0/, XMAX /2.0/, YMIN /-2.0/, YMAX /2.0/, MU /0.3/   MAN  120
      DX = (XMAX-XMIN)/FLOAT(NX-1)                                      MAN  130
      DY = (YMAX-YMIN)/FLOAT(NY-1)                                      MAN  140
      XCUR = 1.0                                                        MAN  150
      YCUR = 1.0                                                        MAN  160
      IF (MOD(NX,2).NE.0) YCUR = FLOAT(NY)                              MAN  170
      IF (MOD(NY,2).NE.0) XCUR = FLOAT(NX)                              MAN  180
      X = XMIN - DX                                                     MAN  190
      DO 20 I=1,NX                                                      MAN  200
        Y = YMIN - DY                                                   MAN  210
        X = X + DX                                                      MAN  220
        DO 10 J=1,NY                                                    MAN  230
          Y = Y + DY                                                    MAN  240
          Z(I,J) = (1.0-MU)*(2.0/SQRT((X-MU)**2+Y**2)+(X-MU)**2+Y**2)   MAN  250
     *      + MU*(2.0/SQRT((X+1.0-MU)**2+Y**2)+(X+1.0-MU)**2+Y**2)      MAN  260
   10   CONTINUE                                                        MAN  270
   20 CONTINUE                                                          MAN  280
      CALL GCONTR(Z, 51, NX, NY, C, NF, 1.E6, WORK, DRAW)               MAN  290
      STOP                                                              MAN  300
      END                                                               MAN  310
      REAL Z(51,51), C(10), CVAL(10), MU                                MAN   10
      INTEGER WORK(1680), L(10), CLAB(10)                               MAN   20
C     DIMENSION OF WORK IS LARGE ENOUGH TO CONTAIN                      MAN   30
C     2*(DIMENSION OF C)*(TOTAL DIMENSION OF Z) USEFUL BITS.  SEE THE   MAN   40
C     BITMAP ROUTINES ACCESSED BY GCONTR.                               MAN   50
      EXTERNAL DRAW                                                     MAN   60
      COMMON /GCTCOM/ XCUR, YCUR, XL, YL, CVAL, CLAB, NCH               MAN   70
      DATA C(1), C(2), C(3), C(4), C(5) /3.05,3.2,3.5,3.50135,3.6/      MAN   80
      DATA C(6), C(7), C(8), C(9), C(10) /3.766413,4.0,4.130149,5.0,    MAN   90
     *  10.0/                                                           MAN  100
      DATA L(1), L(2), L(3), L(4), L(5) /1HA,1HB,1HC,1HD,1HE/           MAN  110
      DATA L(6), L(7), L(8), L(9), L(10) /1HF,1HG,1HH,1HI,1HJ/          MAN  120
      DATA NX /51/, NY /51/, NF /10/, NXG /5/, NYG /5/                  MAN  130
      DATA XMIN /-2.0/, XMAX /2.0/, YMIN /-2.0/, YMAX /2.0/, MU /0.3/   MAN  140
      DATA XLEN /8.0/, YLEN /8.0/                                       MAN  150
C     INITIALIZE PLOTTING SUBROUTINES.                                  MAN  160
      CALL PLOTS                                                        MAN  170
      DX = (XMAX-XMIN)/FLOAT(NX-1)                                      MAN  180
      DY = (YMAX-YMIN)/FLOAT(NY-1)                                      MAN  190
      XL = XLEN/FLOAT(NX)                                               MAN  200
      YL = YLEN/FLOAT(NY)                                               MAN  210
      XCUR = 1.0                                                        MAN  220
      YCUR = 1.0                                                        MAN  230
      IF (MOD(NX,2).NE.0) YCUR = FLOAT(NY)                              MAN  240
      IF (MOD(NY,2).NE.0) XCUR = FLOAT(NX)                              MAN  250
      X = XMIN - DX                                                     MAN  260
      DO 20 I=1,NX                                                      MAN  270
        Y = YMIN - DY                                                   MAN  280
        X = X + DX                                                      MAN  290
        DO 10 J=1,NY                                                    MAN  300
          Y = Y + DY                                                    MAN  310
C     EVALUATE FUNCTION TO BE PLOTTED.                                  MAN  320
          Z(I,J) = (1.0-MU)*(2.0/SQRT((X-MU)**2+Y**2)+(X-MU)**2+Y**2)   MAN  330
     *      + MU*(2.0/SQRT((X+1.0-MU)**2+Y**2)+(X+1.0-MU)**2+Y**2)      MAN  340
   10   CONTINUE                                                        MAN  350
   20 CONTINUE                                                          MAN  360
      DO 30 I=1,NF                                                      MAN  370
        CVAL(I) = C(I)                                                  MAN  380
        CLAB(I) = L(I)                                                  MAN  390
   30 CONTINUE                                                          MAN  400
      NCH = 1                                                           MAN  410
C     PEN UP MOVE TO BELOW LOWER LEFT CORNER OF PAGE.                   MAN  420
C     THIS CALL WORKS DIFFERENTLY ON DIFFERENT MACHINES.  YOU MAY       MAN  430
C     NEED TO CHANGE IT.                                                MAN  440
      CALL PLOT(0.0, -11.0, -3)                                         MAN  450
C     PEN UP MOVE TO 1 INCH ABOVE LOWER LEFT CORNER OF PAGE.            MAN  460
      CALL PLOT(0.0, 1.0, -3)                                           MAN  470
      SX = 8.0/FLOAT(NXG)                                               MAN  480
      SY = 8.0/FLOAT(NXG)                                               MAN  490
C     DRAW A GRID.                                                      MAN  500
      CALL CGRID(1, NXG, SX, 0.0, 0.0, NYG, SY, 0.0, 0.0)               MAN  510
C     DRAW THE CONTOUR PLOTS.                                           MAN  520
      CALL GCONTR(Z, 51, NX, NY, CVAL, NF, 1.0E6, WORK, DRAW)           MAN  530
      XX = 9.0                                                          MAN  540
      YY = 8.0                                                          MAN  550
C     WRITE A TABLE OF CONTOUR LABELS AND VALUES.                       MAN  560
      CALL SYMBOL(XX, YY+0.14, 0.07, 10HCONTOUR ID, 0.0, 10)            MAN  570
      DO 40 I=1,NF                                                      MAN  580
        CALL SYMBOL(XX, YY, 0.07, L(I), 0.0, 2)                         MAN  590
        CALL NUMBER(XX+0.12, YY, 0.07, C(I), 0.0, 5)                    MAN  600
        YY = YY - 0.14                                                  MAN  610
   40 CONTINUE                                                          MAN  620
C     PEN UP MOVE TO BELOW LOWER RIGHT CORNER OF PAGE.                  MAN  630
C     THIS CALL WORKS DIFFERENTLY ON DIFFERENT MACHINES.  YOU MAY NEED  MAN  640
C     TO CHANGE IT, OR YOU MAY NOT NEED IT.                             MAN  650
      CALL PLOT(10.0, -11.0, -3)                                        MAN  660
C     REDUCE PICTURE SIZE, PLOT END OF FILE INFORMATION.                MAN  670
C     THE END OF FILE INFORMATION MAY NOT BE AVAILABLE AT ALL SITES.    MAN  680
C     IF NOT AVAILABLE, CHANGE THE NEXT TWO STATEMENTS TO COMMENTS.     MAN  690
      CALL FACTOR(0.3)                                                  MAN  700
      CALL PLOT(0.0, 0.0, 999)                                          MAN  710
      STOP                                                              MAN  720
C                                                                       MAN  730
      END                                                               MAN  740
      SUBROUTINE DRAW(X, Y, IFLAG)                                      DRA   10
C     THIS SUBROUTINE USES CALCOMP PLOT ROUTINES TO DRAW LINES FOR THE
C     CONTOUR PLOTTING ROUTINE GCONTR.
      REAL CVAL(10)
      INTEGER CLAB(10)
      COMMON /GCTCOM/ XCUR, YCUR, XL, YL, CVAL, CLAB, NCH
      DATA IBLANK /1H /
      IH = IFLAG/10
      IL = IFLAG - 10*IH
      IF (IL.EQ.6) GO TO 40
      IPEN = 2
      IF (IL.EQ.2) IPEN = 3
      IF (IL.EQ.3) IPEN = 3
      XCUR = X
      YCUR = Y
      XX = (X-1.0)*XL
      YY = (Y-1.0)*YL
      CALL PLOT(XX, YY, IPEN)
      IF (IL.LT.2) GO TO 30
      IF (IL.GT.4) GO TO 30
      IF (NCH.LT.1) GO TO 30
      IF (CLAB(IH).EQ.IBLANK) GO TO 30
      IF (CLAB(IH).NE.0) GO TO 10
      CALL NUMBER(XX, YY-0.03, 0.07, CVAL(IH), 0.0, -1)
      GO TO 20
   10 CALL SYMBOL(XX, YY-0.03, 0.07, CLAB(IH), 0.0, NCH)
   20 CALL PLOT(XX, YY, 3)
   30 RETURN
   40 X = XCUR
      Y = YCUR
      RETURN
C
      END
      SUBROUTINE DRAW(X, Y, IFLAG)                                      DRA   10
C
C     DO OUTPUT FOR GCONTR.
C
      INTEGER PRINT
      COMMON /CUR/ XCUR, YCUR
      DATA PRINT /6/
C     PRINT IS THE SYSTEM PRINTER FORTRAN I/O UNIT NUMBER.
      ICONT = IFLAG/10
      JUMP = MOD(IFLAG,10)
      GO TO (10, 20, 30, 40, 50, 60), JUMP
   10 WRITE (PRINT,99999) ICONT, X, Y
      GO TO 70
   20 WRITE (PRINT,99998) ICONT, X, Y
      GO TO 70
   30 WRITE (PRINT,99997) ICONT, X, Y
      GO TO 70
   40 WRITE (PRINT,99996) ICONT, X, Y
      GO TO 70
   50 WRITE (PRINT,99995) ICONT, X, Y
      GO TO 70
   60 WRITE (PRINT,99994)
      X = XCUR
      Y = YCUR
   70 RETURN
99999 FORMAT (17H CONTINUE CONTOUR, I3, 3H TO, 1P2E14.7)
99998 FORMAT (14H START CONTOUR, I3, 19H ON THE BOUNDARY AT, 1P2E14.7)
99997 FORMAT (14H START CONTOUR, I3, 19H IN THE INTERIOR AT, 1P2E14.7)
99996 FORMAT (15H FINISH CONTOUR, I3, 19H ON THE BOUNDARY AT, 1P2E14.7)
99995 FORMAT (15H FINISH CONTOUR, I3, 19H IN THE INTERIOR AT, 1P2E14.7)
99994 FORMAT (33H REQUEST FOR CURRENT PEN POSITION)
      END
      SUBROUTINE CGRID(NOPT, NX, SX, XS, XF, NY, SY, YS, YF)            CGR   10
C
C SUBROUTINE WHICH DRAWS A FRAME AROUND THE PLOT AND DRAWS
C EITHER TICK MARKS OR GRID LINES.
C
C PARAMETERS:  NOPT -- =0, DRAW TICKS ONLY
C                      =1, DRAW GRID LINES
C                      =2, DRAW GRID LINES TO EDGE OF FRAME.
C              NX -- NUMBER OF INTERVALS IN X DIRECTION
C              SX -- SPACING IN INCHES BETWEEN TICK MARKS OR GRID LINES
C                    ALONG THE X AXIS
C              XS -- LOCATION OF FIRST TICK OR GRID LINE ON X AXIS
C              XF -- LOCATION OF RIGHT EDGE OF FRAME
C              NY -- NUMBER OF INTERVALS IN Y DIRECTION
C              SY -- SPACING IN INCHES BETWEEN TICK MARKS OR GRID LINES
C                    ALONG THE Y AXIS
C              YS -- LOCATION OF FIRST TICK OR GRID LINE ON Y AXIS
C              YF -- LOCATION OF TOP EDGE OF FRAME
C ASSUMPTIONS: NX, SX, NY, SY ALL POSITIVE.
C              THE LOWER LEFT-HAND CORNER OF THE FRAME IS DRAWN AT (0,0)
C              IF XS<0, USE 0; IF YS<0, USE 0
C              IF XF<=0, USE NX*SX; IF YF<=0, USE NY*SY.
C
      XINC = SX
      YINC = SY
      XLGTH = FLOAT(NX)*SX
      YLGTH = FLOAT(NY)*SY
      XMIN = AMAX1(XS,0.0)
      YMIN = AMAX1(YS,0.0)
      XMAX = AMAX1(XF,XLGTH+XMIN)
      YMAX = AMAX1(YF,YLGTH+YMIN)
C
C     DRAW FRAME.
C
      CALL PLOT(0.0, 0.0, 3)
      CALL PLOT(XMAX, 0.0, 2)
      CALL PLOT(XMAX, YMAX, 2)
      CALL PLOT(0.0, YMAX, 2)
      CALL PLOT(0.0, 0.0, 2)
      IF (NOPT.NE.0) GO TO 130
C
C     DRAW TICK MARKS.
C
      DO 120 J=1,4
        GO TO (10, 50, 20, 40), J
   10   X2 = 0.0
        IF (XMIN.NE.0.0) X2 = XMIN - SX
        Y2 = 0.0
        GO TO 30
   20   XINC = -SX
        X2 = XMIN + XLGTH + SX
        IF (XMAX.EQ.XMIN+XLGTH) X2 = XMAX
        Y2 = YMAX
   30   Y1 = Y2
        Y2 = Y2 + SIGN(0.125,XINC)
        N = NX
        IF (ABS(XMAX-XMIN-XLGTH)+ABS(XMIN)) 70, 80, 70
   40   YINC = -SY
        Y2 = YMIN + YLGTH + SY
        IF (YMAX.EQ.YMIN+YLGTH) Y2 = YMAX
        X2 = 0.0
        GO TO 60
   50   Y2 = 0.0
        IF (YMIN.NE.0.0) Y2 = YMIN - SY
        X2 = XMAX
   60   X1 = X2
        N = NY
        X2 = X2 - SIGN(0.125,YINC)
        IF (ABS(YMAX-YMIN-YLGTH)+ABS(YMIN)) 70, 80, 70
   70   N = N + 1
   80   DO 110 I=1,N
          IF (MOD(J,2).EQ.0) GO TO 90
          X2 = X2 + XINC
          X1 = X2
          GO TO 100
   90     Y2 = Y2 + YINC
          Y1 = Y2
  100     CALL PLOT(X1, Y1, 3)
          CALL PLOT(X2, Y2, 2)
  110   CONTINUE
  120 CONTINUE
      GO TO 240
C
C     DRAW GRID LINES
C
  130 X1 = XMIN
      X2 = XMIN + XLGTH
      IF (NOPT.NE.2) GO TO 140
      X1 = 0.0
      X2 = XMAX
  140 Y1 = YMIN - SY
      N = NY + 1
      IF (YMAX.EQ.YMIN+YLGTH) N = N - 1
      IF (YMIN.NE.0.0) GO TO 150
      Y1 = 0.0
      N = N - 1
  150 IF (N.LE.0) GO TO 170
      J = 1
      DO 160 I=1,N
        J = -J
        Y1 = Y1 + SY
        CALL PLOT(X1, Y1, 3)
        CALL PLOT(X2, Y1, 2)
        XX = X1
        X1 = X2
        X2 = XX
  160 CONTINUE
  170 Y1 = YMIN + YLGTH
      Y2 = YMIN
      IF (NOPT.NE.2) GO TO 180
      Y1 = YMAX
      Y2 = 0.0
  180 N = NX + 1
      IF (J.LT.0) GO TO 200
      X1 = XMIN - SX
      IF (XMAX.EQ.XMIN+XLGTH) N = N - 1
      IF (XMIN.NE.0.0) GO TO 190
      X1 = 0.0
      N = N - 1
  190 IF (N.LE.0) GO TO 240
      XINC = SX
      GO TO 220
  200 X1 = XMIN + XLGTH + SX
      IF (XMIN.EQ.0.0) N = N - 1
      IF (XMAX.NE.XLGTH+XMIN) GO TO 210
      N = N - 1
      X1 = XMAX
  210 XINC = -SX
  220 DO 230 I=1,N
        X1 = X1 + XINC
        CALL PLOT(X1, Y1, 3)
        CALL PLOT(X1, Y2, 2)
        XX = Y1
        Y1 = Y2
        Y2 = XX
  230 CONTINUE
  240 RETURN
C
      END
