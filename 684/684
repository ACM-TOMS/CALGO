C      ALGORITHM 684, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 16, NO. 3, PP. 253-257.
      PROGRAM  DFC1C2
C
C     EXAMPLE FOR COMPUTING VALUES Z(I,J)
C     FOR RECTANGULAR GRIDS XX(I), YY(J), I=1,NXI, J=1,NYI
C     FROM SCATTERED DATA XD(K),YD(K),ZD(K), K=1,ND
C     USING A C2 NONIC AND C1 QUINTIC
C     POLYNOMIAL REPRESENTATION OVER TRIANGLES.
C
C     CONTOUR PLOTS ARE PRESENTED IN THE CORRESPONDING PAPER
C
C     AUTHOR      : A. PREUSSER, 1989
C                   FRITZ-HABER-INSTITUT
C                   DER MAX-PLANCK-GESELLSCHAFT
C                   FARADAYWEG 4-6
C                   D-1000 BERLIN 33
C
C
*D    DOUBLE PRECISION Z
      DIMENSION  XD(30),YD(30),ZD(30)
      COMMON XX(101),YY(081),Z(101,081),Z1(101,081)
      DIMENSION IWK(210),WK(30,14)
      DATA  ND/30/
      DATA  NXI,NYI /101,081/
C     WORKSPACE LENGTH IWK ND*7, WK ND*14
C
C     XD,YD,ZD     - COORDINATES OF ND DATA POINTS
C     XX,YY        - COORDINATES OF A NXI*NYI RECTANGULAR GRID
C     Z,Z1         - ARRAYS CONTAINING THE INTERPOLATED VALUES
C
C
C     THE XD- AND YD-VALUES ARE
C     TAKEN FROM ACM ALGORITHM 526 BY H.AKIMA
C
      DATA  XD(1), XD(2), XD(3), XD(4), XD(5), XD(6),
     1      XD(7), XD(8), XD(9), XD(10),XD(11),XD(12),
     2      XD(13),XD(14),XD(15),XD(16),XD(17),XD(18),
     3      XD(19),XD(20),XD(21),XD(22),XD(23),XD(24),
     4      XD(25),XD(26),XD(27),XD(28),XD(29),XD(30)/
     5      11.16, 24.20, 19.85, 10.35, 19.72,  0.00,
     6      20.87, 19.99, 10.28,  4.51,  0.00, 16.70,
     7       6.08, 25.00, 14.90,  0.00,  9.66,  5.22,
     8      11.77, 15.10, 25.00, 25.00, 14.59, 15.20,
     9       5.23,  2.14,  0.51, 25.00, 21.67,  3.31/
      DATA  YD(1), YD(2), YD(3), YD(4), YD(5), YD(6),
     1      YD(7), YD(8), YD(9), YD(10),YD(11),YD(12),
     2      YD(13),YD(14),YD(15),YD(16),YD(17),YD(18),
     3      YD(19),YD(20),YD(21),YD(22),YD(23),YD(24),
     4      YD(25),YD(26),YD(27),YD(28),YD(29),YD(30)/
     5       1.24, 16.23, 10.72,  4.11,  1.39, 20.00,
     6      20.00,  4.62, 15.16, 20.00,  4.48, 19.65,
     7       4.58, 11.87,  3.12,  0.00, 20.00, 14.66,
     8      10.47, 17.19,  3.87,  0.00,  8.71,  0.00,
     9      10.72, 15.03,  8.37, 20.00, 14.36,  0.13/
C
C
C     SET GRID COORDINATES
      DX= 25./(NXI-1)
      DY= 20./(NYI-1)
      DO 10 I=1,NXI
   10   XX(I)= (I-1)*DX
      DO 20 I=1,NYI
   20   YY(I)= (I-1)*DY
C
C     SET ZD-VALUES
      DO 30 I=1,ND
   30   ZD(I)= 5 + MOD(I,2)*10
C
C     COMPUTE RECTANGULAR GRIDS
      CALL C2GRID (ND,XD,YD,ZD,NXI,NYI,XX,YY,Z,NXI,WK,14*ND,
     1             IWK,7*ND)
C
      CALL C1GRID (ND,XD,YD,ZD,NXI,NYI,XX,YY,Z1,NXI,WK,6*ND,
     1             IWK,7*ND)
C
C     COMPUTE DIFFERENCE
      DO 40 J=1,NYI
        DO 40 I=1,NXI
          Z(I,J)= Z(I,J) - Z1(I,J)
  40  CONTINUE
C
C     THE FOLLOWING STATEMENT MAY BE ACTIVATED FOR A CONTOUR PLOT.
C     FOR A DESCRIPTION OF SUBROUTINE FARBE,
C     SEE ACM ALGORITHM 671, ACM TOMS 15, MARCH 1989, BY A.PREUSSER.
C     THE ACCOMPANYING PAPER TO THE C1-C2 ROUTINES ALSO CONTAINS
C     A CONTOUR PLOT OF THE DIFFERENCE SURFACE.
*     CALL FARBE(Z,NXI,NXI,NYI,2)
C
C     THE FOLLOWING WRITE STATEMENTS CAN BE TAKEN AS A ROUGH
C     CHECK FOR CORRECT INSTALLATION.
C     NORMALLY, UP TO THREE DIGITS CAN BE LOST BECAUSE OF
C     ROUNDING ERRORS.
C     ROUNDING ERRORS WILL BE REDUCED BY ACTIVATING
C     THE STATEMENTS WITH   *D   IN THE FIRST TWO COLUMNS.
C     THEN DOUBLE PRECISION WILL BE USED FOR THE
C     COMPUTATIONS OF THE C2-POLYNOMIALS.
C     SINCE A HIGH PRECISION COMPUTATION ONLY MAKES SENSE
C     WITH HIGH PRECISION DATA, THE PARAMETERS E AND L
C     (EPS AND NIT) FOR THE COMPUTATION OF THE PARTIAL
C     DERIVATIVES BY SUBROUTINE GRANKA SHOULD ALSO BE ADJUSTED
C     IN THIS CASE.
C     (L=15 IS NORMALLY SUFFICIENT FOR 7 DIGITS).
      WRITE (*,999) 3.035646151154D0,Z(44,03)
  999 FORMAT (' EXPECTED VALUE',F20.12/
     .        ' COMPUTED VALUE',F20.12)
C
C
C     ALTERNATIVE WAY
C     FOR GETTING THE SAME RESULTS
C
C     INITIALIZATION FOR C2PNT
      CALL C2INI (ND,XD,YD,ZD,WK,14*ND,IWK,7*ND)
C
C     LOOP FOR COMPUTATION OF POINTS WITH C2PNT
      IST= 1
      DO 120 J=1,NYI
        DO 110 I=1,NXI
          CALL C2PNT (ND,XD,YD,ZD,XX(I),YY(J),WK,IWK,
     1                 IST,Z(I,J),IOUTS)
  110   CONTINUE
  120 CONTINUE
C
C     INITIALIZATION FOR C1PNT
      CALL C1INI (ND,XD,YD,ZD,WK,6*ND,IWK,7*ND)
C
C     LOOP FOR COMPUTATION OF POINTS WITH C1PNT
      IST= 1
      DO 160 J=1,NYI
        DO 150 I=1,NXI
          CALL C1PNT (ND,XD,YD,ZD,XX(I),YY(J),WK,IWK,
     1                 IST,Z1(I,J),IOUTS)
  150   CONTINUE
  160 CONTINUE
C
C     DIFFERENCE
      DO 220 J=1,NYI
        DO 210 I=1,NXI
          Z(I,J)= Z(I,J) - Z1(I,J)
  210   CONTINUE
  220 CONTINUE
C
      WRITE (*,999) 3.035646151154D0,Z(44,03)
C
      STOP
      END
      SUBROUTINE C2GRID (ND,XD,YD,ZD,NX,NY,XX,YY,Z,NDIMX,WK,NW,
     1                   IWK,IW)
C
*D    DOUBLE PRECISION Z
      DIMENSION XD(ND),YD(ND),ZD(ND)
      DIMENSION XX(NX),YY(NY)
      DIMENSION Z(NDIMX,NY)
      DIMENSION WK(ND,14),IWK(IW)
C
C**********************************************************************
C
C     COMPUTATION OF A RECTANGULAR GRID FROM SCATTERED DATA
C     =========================================================
C     USING A C2 NONIC POLYNOMIAL REPRESENTATION OVER TRIANGLES
C
C
C     AUTHOR      : A. PREUSSER, 1989
C                   FRITZ-HABER-INSTITUT
C                   DER MAX-PLANCK-GESELLSCHAFT
C                   FARADAYWEG 4-6
C                   D-1000 BERLIN 33
C
C     INPUT PARAMETERS
C     ----------------
C
C     ND          NUMBER OF DATA POINTS
C     XD,YD,ZD    ARRAYS OF LENGTH ND  WITH
C                 X,Y,Z-COORDINATES OF DATA POINTS
C     NX,NY       NUMBER OF GRID LINES IN X- AND
C                 Y-DIRECTION.
C     XX,YY       ARRAYS OF LENGTH NX, RESP. NY,
C                 WITH COORDINATES FOR GRID LINES
C                 X= XX(I), Y= YY(J), I=1,NX; J=1,NY,
C                 IN ASCENDING ORDER.
C     Z           OUTPUT PARAMETER (SEE BELOW).
C     NDIMX       FIRST DIMENSION OF THE TWO-DIMENSIONAL
C                 ARRAY Z AS DECLARED IN THE CALLING
C                 (SUB-)PROGRAM.
C     WK          REAL WORK ARRAY OF LENGTH NW.
C     NW          DIMENSION OF ARRAY WK,
C                 NW.GE.14*ND.
C     IWK         INTEGER WORK ARRAY OF LENGTH IW.
C     IW          DIMENSION OF ARRAY IWK,
C                 IW.GE.7*ND.
C
C
C     OUTPUT PARAMETER
C     ----------------
C
C     Z           TWO-DIMENSIONAL ARRAY WITH DIMENSION
C                 (NDIMX,NY) CONTAINING THE INTERPOLATED
C                 Z-VALUES AT THE INTERSECTIONS OF THE
C                 GRID LINES (=GRID POINTS).
C                 Z(I,J) CORRESPONDS TO THE POINT WITH
C                 THE COORDINATES XX(I),YY(J).
C                 IF XX(I),YY(J) IS OUTSIDE THE TRIANGULATION
C                 OF THE ND DATA POINTS, Z(I,J) IS NOT ALTERED.
C
C     METHODS
C     -------
C
C     A CONVEX TRIANGULAR MESH IS FORMED IN THE X-Y PLANE,
C     WITH THE PROJECTIONS OF THE ND DATA POINTS AS NODES.
C     (SUBROUTINE TRMESH FROM ACM ALG. 624 BY R.RENKA)
C
C     PARTIAL DERIVATIVES UP TO FOURTH ORDER ARE ESTIMATED
C     AT THE NODES.
C     (SUBSEQUENT CALLS TO GRANKA, A MODIFIED VERSION OF GRADG
C      FROM ACM ALG. 624 BY R.RENKA)
C
C     A BIVARIATE NONIC HERMITE POLYNOMIAL IS FORMED FOR
C     EACH TRIANGLE.
C
C     THE POLYNOMIAL IS EVALUATED BY HORNER'S METHOD, IN CASE
C     THE REQUESTED POINT XX(I),YY(J) IS INSIDE THE TRIANGULATION.
C     IF THE POINT IS OUTSIDE, Z(I,J) IS NOT CHANGED.
C
C***********************************************************************
C
      DIMENSION XL(3),YL(3),ZL(3),PDL(14,3)
C
C     ERROR CHECKS
C
      IF (ND.LT.3) THEN
        WRITE (*,901) ND
        GOTO 1000
       END IF
      IF (NX.LT.1) THEN
        WRITE (*,902) NX
        GOTO 1000
       END IF
      IF (NY.LT.1) THEN
        WRITE (*,903) NY
        GOTO 1000
       END IF
      IF (NX.GT.NDIMX) THEN
        WRITE (*,904) NX,NDIMX
        GOTO 1000
       END IF
      IF (NW.LT.14*ND)  THEN
        WRITE (*,905) NW,14*ND
        GOTO 1000
       END IF
      IF (IW.LT.7*ND)  THEN
        WRITE (*,906) IW,7*ND
        GOTO 1000
       END IF
      IXERR= 0
      DO 20 I=2,NX
        IF (XX(I)-XX(I-1).GT.0.) GOTO 20
        IXERR= 1
   20 CONTINUE
      IF (IXERR.EQ.1)  WRITE (*,907) I,XX(I),I-1,XX(I-1)
      IYERR= 0
      DO 30 I=2,NY
        IF (YY(I)-YY(I-1).GT.0.) GOTO 30
        IYERR= 1
   30 CONTINUE
      IF (IYERR.EQ.1)  WRITE (*,908) I,YY(I),I-1,YY(I-1)
  901 FORMAT (' ***ERROR*** ND.LT.3, ND=',I10)
  902 FORMAT (' ***ERROR*** NX.LT.1, NX=',I10)
  903 FORMAT (' ***ERROR*** NY.LT.1, NY=',I10)
  904 FORMAT (' ***ERROR*** NX.GT.NXDIM. NX,NXDIM=',I10,',',I10)
  905 FORMAT (' ***ERROR*** NW.LT.14*ND. NW,14*ND=',I10,',',I10)
  906 FORMAT (' ***ERROR*** IW.LT.7*ND. IW,7*ND=',I10,',',I10)
  907 FORMAT ('0***ERROR***          '/
     A        '      X-COORDINATE', I5,'=',E15.7/
     B        ' .LE. X-COORDINATE', I5,'=',E15.7)
  908 FORMAT ('0***ERROR***          '/
     A        '      Y-COORDINATE', I5,'=',E15.7/
     B        ' .LE. Y-COORDINATE', I5,'=',E15.7)
C
C     TRIANGULATION AND COMPUTATION OF PARTIALS
      CALL C2INI (ND,XD,YD,ZD,WK,NW,IWK,IW)
*
C     COMPUTE INTERPOLATED VALUES  Z(I,J)
C
C     LOOPS TRIANGLES
      INDF= 1
      N21= 6*ND-1
      DO 500 K1=1,ND-2
        I1= K1
        INDL= IWK(N21+K1)
        DO 400 INDX=INDF,INDL
          I2= IWK(INDX)
          I3= IWK(INDX+1)
          IF (INDX.EQ.INDL) I3= IWK(INDF)
          IF (I2.LT.K1 .OR. I3.LT.K1) GOTO 400
C
C         START COMPUTATIONS FOR TRIANGLE
C
C         EXTREME VALUES OF VERTEX COORDINATES
          XMIN= AMIN1(XD(I1),XD(I2),XD(I3))
          XMAX= AMAX1(XD(I1),XD(I2),XD(I3))
          YMIN= AMIN1(YD(I1),YD(I2),YD(I3))
          YMAX= AMAX1(YD(I1),YD(I2),YD(I3))
C
C         FIND MIN AND MAX INDICES IN X OF GRID
C
          IXMIN= 0
          DO 110 I=1,NX
            IF (XX(I).LT.XMIN) GOTO 110
            IXMIN= I
            GOTO 115
  110     CONTINUE
          IF (IXMIN.EQ.0) GOTO 400
  115     CONTINUE
C
          IXMAX= 0
          DO 120 I=NX,1,-1
            IF (XX(I).GT.XMAX) GOTO 120
            IXMAX= I
            GOTO 125
  120     CONTINUE
          IF (IXMAX.EQ.0) GOTO 400
  125     CONTINUE
C
C
C         LOOPS GRID POINTS
          ITFLAG= 0
          DO 390 J=1,NY
            IF (YY(J).LT.YMIN .OR. YY(J).GT.YMAX) GOTO 390
            DO 370 I=IXMIN,IXMAX
C
C             START COMPUTATIONS FOR GRID POINT
C
C             CHECK IF GRID POINT INSIDE TRIANGLE
              IF ((XD(I1)-XX(I))*(YD(I2)-YY(J)) -
     1            (YD(I1)-YY(J))*(XD(I2)-XX(I)) .LT. 0. .OR.
     2            (XD(I2)-XX(I))*(YD(I3)-YY(J)) -
     3            (YD(I2)-YY(J))*(XD(I3)-XX(I)) .LT. 0. .OR.
     4            (XD(I3)-XX(I))*(YD(I1)-YY(J)) -
     5            (YD(I3)-YY(J))*(XD(I1)-XX(I)) .LT.0) GOTO 370
C
C             CHECK IF COEFFICIENTS FOR TRIANGLE ALREADY COMPUTED
              IF (ITFLAG.NE.0) GOTO 350
              ITFLAG= 1
C
C             LOAD COORDINATES
              XL(1)= XD(I1)
              YL(1)= YD(I1)
              XL(2)= XD(I2)
              YL(2)= YD(I2)
              XL(3)= XD(I3)
              YL(3)= YD(I3)
C
C             IF THE FOLLOWING 3 GROUPS OF STATEMENTS WITH
C             AN ASTERIX IN THE FIRST COLUMN ARE ACTIVATED,
C             RESULTS DO NOT CHANGE, ROUNDING ERRORS EXCEPTED.
*             XL(2)= XD(I1)
*             YL(2)= YD(I1)
*             XL(3)= XD(I2)
*             YL(3)= YD(I2)
*             XL(1)= XD(I3)
*             YL(1)= YD(I3)
C
C             LOAD COORDINATES
              ZL(1)= ZD(I1)
              ZL(2)= ZD(I2)
              ZL(3)= ZD(I3)
*             ZL(2)= ZD(I1)
*             ZL(3)= ZD(I2)
*             ZL(1)= ZD(I3)
*
C             LOAD PARTIAL DERIVATIVES
              DO 300 K=1,14
                PDL(K,1)= WK(I1,K)
                PDL(K,2)= WK(I2,K)
                PDL(K,3)= WK(I3,K)
*               PDL(K,2)= WK(I1,K)
*               PDL(K,3)= WK(I2,K)
*               PDL(K,1)= WK(I3,K)
  300         CONTINUE
C
C             COMPUTE COEFFICIENTS ALONG SIDES
              CALL C2SIDE (XL,YL,ZL,PDL,2)
C
C             COMPUTE COEFFICIENTS FOR INSIDE OF TRIANGLE
C             AND CONSTANTS FOR TRANSF. CARTESIAN TO TRIANGULAR
              CALL C2INSD (PDL)
C
C             EVALUATE POLYNOMIAL AT XP,YP
  350         CONTINUE
              XR= XX(I) -XL(3)
              YR= YY(J) -YL(3)
              CALL C2HORN (XR,YR,Z(I,J))
C
  370       CONTINUE
  390     CONTINUE
  400   CONTINUE
        INDF= INDL + 1
  500 CONTINUE
C
 1000 RETURN
      END
      SUBROUTINE C2INI (ND,XD,YD,ZD,WK,NW,IREN,IR)
C
      DIMENSION XD(ND),YD(ND),ZD(ND)
      DIMENSION WK(ND,14),IREN(IR)
C
C**********************************************************************
C
C     TRIANGULATE DATA AND COMPUTE PARTIAL DERIVATIVES UP TO 4TH ORDER
C     ================================================================
C
C
C     AUTHOR      : A. PREUSSER, 1989
C                   FRITZ-HABER-INSTITUT
C                   DER MAX-PLANCK-GESELLSCHAFT
C                   FARADAYWEG 4-6
C                   D-1000 BERLIN 33
C
C     INPUT PARAMETERS
C     ----------------
C
C     ND          NUMBER OF DATA POINTS
C     XD,YD,ZD    ARRAYS OF LENGTH ND  WITH
C                 X,Y,Z-COORDINATES OF DATA POINTS
C
C     OUTPUT PARAMETERS
C     -----------------
C
C     WK          REAL ARRAY OF LENGTH NW
C                 CONTAINING THE PARTIAL DERIVATIVES.
C                 FIRST INDEX FOR POINT, SECOND FOR DERIVATIVE.
C                 INDICES FOR DERIVATIVES:
C
C                              1    2
C                             ZX    ZY
C                          3     4      5
C                        ZXX    ZXY     ZYY
C                      6      7      8      9
C                     ZXXX   ZXXY   ZXYY   ZYYY
C                   10     11     12     13     14
C                  ZXXXX  ZXXXY  ZXXYY  ZXYYY  ZYYYY
C
C
C     NW          DIMENSION OF ARRAY WK,
C                 NW.GE.14*ND.
C     IREN        INTEGER ARRAY OF LENGTH IR
C                 CONTAINING THE DATA STRUCTURE FOR THE
C                 DESCRIPTION OF THE TRIANGELS AS
C                 DEFINED BY R.RENKA FOR ACM ALG. 624
C                 (ARRAYS IADJ AND IEND).
C     IR          DIMENSION OF ARRAY IREN,
C                 IR.GE.7*ND.
C
C
C     METHODS
C     -------
C
C     A CONVEX TRIANGULAR MESH IS FORMED IN THE X-Y PLANE,
C     WITH THE PROJECTIONS OF THE ND DATA POINTS AS NODES.
C     (SUBROUTINE TRMESH FROM ACM ALG. 624 BY R.RENKA)
C
C     PARTIAL DERIVATIVES UP TO FOURTH ORDER ARE ESTIMATED
C     AT THE NODES.
C     (SUBSEQUENT CALLS TO GRANKA, A MODIFIED VERSION OF GRADG
C      FROM ACM ALG. 624 BY R.RENKA)
C
C***********************************************************************
C
C     ERROR CHECKS
C
      IF (ND.LT.3) THEN
        WRITE (*,901) ND
        GOTO 1000
       END IF
      IF (NW.LT.14*ND)  THEN
        WRITE (*,905) NW,14*ND
        GOTO 1000
       END IF
      IF (IR.LT.7*ND)  THEN
        WRITE (*,906) IR,7*ND
        GOTO 1000
       END IF
  901 FORMAT (' ***ERROR*** ND.LT.3, ND=',I10)
  905 FORMAT (' ***ERROR*** NW.LT.14*ND. NW,14*ND=',I10,',',I10)
  906 FORMAT (' ***ERROR*** IR.LT.7*ND. IR,7*ND=',I10,',',I10)
C
C     TRIANGULATION
      N2= 6*ND
      CALL TRMESH (ND,XD,YD,IREN,IREN(N2),IER)
C     TRMESH IS PART OF ACM ALGORITHM 624 BY R.RENKA
      IF (IER.EQ.3) WRITE (*,911)
      IF (IER.EQ.3) GOTO 1000
  911 FORMAT (' ***ERROR*** ALL POINTS COLLINEAR')
C
C     CHECK FOR DUPLICATE POINTS
      INDF= 1
      N21= N2-1
      DO 100 K1=1,ND-1
        IND1= K1
        INDL= IREN(N21+K1)
        DO 40 INDX=INDF,INDL
          IND2= IREN(INDX)
          IF (IND2.LT.K1) GOTO 40
C         THE FOLLOWING TWO STATEMENTS MAY BE ACTIVATED
C         FOR PLOTTING THE MESH LINES
*         CALL PLOT (XD(IND1),YD(IND1),3)
*         CALL PLOT (XD(IND2),YD(IND2),2)
          IF (ABS(XD(IND1)-XD(IND2))+
     .        ABS(YD(IND1)-YD(IND2)).NE.0.) GOTO 40
          WRITE (*,912) IND1,IND2
          GOTO 1000
   40   CONTINUE
        INDF= INDL + 1
  100 CONTINUE
  912 FORMAT (' ***ERROR*** DUPLICATE POINTS',2I10)
C
C
C     COMPUTE PARTIAL DERIVATIVES BY RENKA'S GLOBAL METHOD
C
C     GRANKA IS A MODIFIED VERSION OF GRADG FROM
C     ACM ALGORITHM 624
C
      DO 140 J=1,14
        DO 140 I=1,ND
  140     WK(I,J)= 0.
C
C     COMPUTE ZX AND ZY
C     L= NUMBER OF ITERATIONS
      L= 15
      E= 0.
      CALL GRANKA (ND,XD,YD,ZD,IREN,IREN(N2),E,L,WK,IER)
C
      L= 15
C     COMPUTE ZXX AN ZXY
      CALL GRANKA (ND,XD,YD,WK,IREN,IREN(N2),E,L,WK(1,3),IER)
C
C     COMPUTE ZYX AND ZYY
      L= 15
      CALL GRANKA (ND,XD,YD,WK(1,2),IREN,IREN(N2),E,L,WK(1,5),IER)
C
C     COMPUTE ZXY= (ZXY+ZYX)/2
      DO 200 I=1,ND
        WK(I,4)= (WK(I,4) + WK(I,5))*0.5
        WK(I,5)= WK(I,6)
        WK(I,6)= 0.
  200 CONTINUE
C
C     COMPUTE ZXXX AND ZXXY
      L= 15
      CALL GRANKA (ND,XD,YD,WK(1,3),IREN,IREN(N2),E,L,WK(1,6),IER)
C
C     COMPUTE ZYYX UND ZYYY
      L= 15
      CALL GRANKA (ND,XD,YD,WK(1,5),IREN,IREN(N2),E,L,WK(1,8),IER)
C
C     COMPUTE ZXXYY ON WK(I,11)
      L= 15
      CALL GRANKA (ND,XD,YD,WK(1,7),IREN,IREN(N2),E,L,WK(1,10),IER)
C
C     COMPUTE ZYYXX ON WK(I,12)
      L= 15
      CALL GRANKA (ND,XD,YD,WK(1,8),IREN,IREN(N2),E,L,WK(1,12),IER)
C
C     COMPUTE ZXXYY= (ZXXYY + ZYYXX)*0.5
      DO 250 I=1,ND
        WK(I,12)= (WK(I,11) + WK(I,12))*0.5
        WK(I,10)= 0.
        WK(I,11)= 0.
        WK(I,13)= 0.
  250 CONTINUE
C
C     COMPUTE ZXXXX AND ZXXXY ON WK(I,10) AND WK(I,11)
      L= 15
      CALL GRANKA (ND,XD,YD,WK(1,6),IREN,IREN(N2),E,L,WK(1,10),IER)
C
C     COMPUTE ZYYYX AND ZYYYY ON WK(I,13) AND WK(I,14)
      L= 15
      CALL GRANKA (ND,XD,YD,WK(1,9),IREN,IREN(N2),E,L,WK(1,13),IER)
C
C
 1000 CONTINUE
      RETURN
      END
      SUBROUTINE C2PNT (ND,XD,YD,ZD,XP,YP,WK,
     .                  IREN,IST,ZP,IOUTS)
*D    DOUBLE PRECISION ZP
      DIMENSION IREN(7*ND)
      DIMENSION XD(ND), YD(ND), ZD(ND), WK(ND,14)
     A,         X(3),Y(3),Z(3),PDL(14,3)
C
C***********************************************************
C
C     INTERPOLATE Z-VALUE OF A POINT IN A C2-REPRESENTATION
C     -----------------------------------------------------
C                   ON A SET OF TRIANGLES
C                   ---------------------
C
C
C     AUTHOR:
C     A.PREUSSER
C     FRITZ-HABER-INSTITUT DER MPG
C     FARADAYWEG 4-6
C     D-1000 BERLIN 33
C
C
C
C
C GIVEN A TRIANGULATION OF A SET OF POINTS IN THE PLANE,
C IN A DATA STRUCTURE AS DEVELOPED BY R. RENKA FOR ACM ALG. 624,
C AND VALUES AND PARTIAL DERIVATIVES AT THESE POINTS,
C THIS ROUTINE COMPUTES A VALUE ZP AT XP,YP.
C
C THE PREREQUISITE VALUES CAN BE COMPUTED BY SUBROUTINE C2INI.
C
C INPUT PARAMETERS
C ----------------
C       ND -    NUMBER OF NODES IN THE MESH.
C               ND .GE. 3.
C               THE VALUE OF ND MUST BE EXACTLY
C               THE SAME AS SPECIFIED FOR THE CALL
C               OF SUBROUTINE C2INI.
C
C       XD,YD - VECTORS OF COORDINATES OF THE
C               NODES IN THE MESH.
C
C          ZD - VECTOR OF DATA VALUES AT THE
C               NODES.
C
C       XP,YP - COORDINATES OF A POINT AT WHICH
C               FUNCTION VALUE ZP IS REQUIRED
C
C          WK - ND BY 14 ARRAY CONTAINING
C               PARTIAL DERIVATIVES AT THE
C               NODES.
C               FIRST INDEX INDICATES POINT,
C               SECOND INDEX DERIVATIVE.
C               ORDER OF COLUMNS:
C                             ZX,ZY
C                         ZXX, ZXY, ZYY
C                    ZXXX, ZXXY, ZXYY, ZYYY
C               ZXXXX, ZXXXY, ZXXYY, ZXYYY, ZYYYY
C
C        IREN - INTEGER ARRAY OF LENGTH 7*ND.
C               CAN BE COMPUTED BY A CALL TO SUBROUTINE
C               C2INI. IT CONTAINS THE DATA FOR
C               THE DEFINITION OF THE TRIANGLES.
C               (ARRAYS IADJ AND IEND OF ALG. 624)
C
C         IST - INDEX OF THE STARTING NODE IN
C               THE SEARCH FOR A TRIANGLE CON-
C               TAINING (XP,YP).  1 .LE. IST
C               .LE. ND.  THE OUTPUT VALUE OF
C               IST FROM A PREVIOUS CALL MAY
C               BE A GOOD CHOICE.
C
C   INPUT PARAMETERS OTHER THAN IST ARE NOT ALTERED BY THIS
C   ROUTINE.
C
C OUTPUT PARAMETERS
C -----------------
C
C         IST - INDEX OF ONE OF THE VERTICES OF
C               THE TRIANGLE CONTAINING (XP,YP),
C
C          ZP - VALUE TO BE COMPUTED AT (XP,YP)
C               =0, IF XP,YP IS OUTSIDE THE TRIANGULATION
C
C       IOUTS - =0, IF (XP,YP) INSIDE THE TRIANGULATION
C                   AND ZP CONTAINS THE INTERPOLATED VALUE.
C               =1, IF (XP,YP) IS OUTSIDE THE TRIANGULATION
C                   AND ZP=0.
C
C     REMARK
C     ------
C     THIS ROUTINE SAVES THE INDICES OF THE THREE VERTICES OF
C     THE TRIANGLE IN WHICH THE POINT XP,YP IS FOUND.
C     IF, IN THE NEXT CALL, THE SAME INDICES ARE FOUND,
C     IT ASSUMES THAT /C2PCO/ AND /C2CCO/ HAVE NOT BEEN
C     ALTERED BY THE USER. THEN THE COEFFICIENTS FOR
C     THE TRIANGLE ARE NOT CALCULATED BUT TAKEN FROM /C2PCO/
C     AND /C2CCO/.
C***********************************************************
C
      SAVE I1OLD,I2OLD,I3OLD,X,Y
      DATA I1OLD,I2OLD,I3OLD /0,0,0/
C
C     ERROR CHECKS
      IF (ND.LT.3) THEN
         WRITE (*,901) ND
         GOTO 3000
       END IF
      IF (IST.LT.1) THEN
         WRITE (*,902) IST
         GOTO 3000
       END IF
      IF (IST.GT.ND) THEN
         WRITE (*,903) IST,ND
         GOTO 3000
       END IF
  901 FORMAT (' ***ERROR*** ND.LT.3, ND= ',I10)
  902 FORMAT (' ***ERROR*** IST.LT.1, IST=',I10)
  903 FORMAT (' ***ERROR*** IST.GT.ND, IST,ND=',I10,',',I10)
C
C     FIND A TRIANGLE CONTAINING P
C
      N2= 6*ND
      CALL TRFIND(IST,XP,YP,XD,YD,IREN,IREN(N2),I1,I2,I3)
      IF (I1 .EQ. 0 .OR. I2.EQ.0 .OR. I3.EQ.0) GOTO 1000
      IST = I1
C
C     CHECK IF SAME TRIANGLE AS BEFORE
C     (COUNTER-CLOCKWISE ORDER OF VERTICES)
      IF (I1.EQ.I1OLD .AND. I2.EQ.I2OLD .AND. I3.EQ.I3OLD
     1.OR.I1.EQ.I2OLD .AND. I2.EQ.I3OLD .AND. I3.EQ.I1OLD
     2.OR.I1.EQ.I3OLD .AND. I2.EQ.I1OLD .AND. I3.EQ.I2OLD)
     3                                           GOTO 500
C
C     LOAD PARTIAL DERIVATIVES
C
      K1= 1
      K2= 2
      K3= 3
C     IF ACTIVATING THE FOLLOWING STATEMENTS,
C     RESULTS SHOULD NOT CHANGE.
*     K1= 2
*     K2= 3
*     K3= 1
      DO 100 J=1,14
        PDL(J,K1)= WK(I1,J)
        PDL(J,K2)= WK(I2,J)
        PDL(J,K3)= WK(I3,J)
  100 CONTINUE
C
C     LOAD COORDINATES
C
      X(K1)= XD(I1)
      Y(K1)= YD(I1)
      X(K2)= XD(I2)
      Y(K2)= YD(I2)
      X(K3)= XD(I3)
      Y(K3)= YD(I3)
      Z(K1)= ZD(I1)
      Z(K2)= ZD(I2)
      Z(K3)= ZD(I3)
C
C     COMPUTE COEFFICIENTS ALONG SIDES
      CALL C2SIDE (X,Y,Z,PDL,2)
C
C     COMPUTE COEFFICIENTS FOR INSIDE OF TRIANGLE
C     AND CONSTANTS FOR TRANSF. CARTESIAN TO TRIANGULAR
      CALL C2INSD (PDL)
C
C     EVALUATE POLYNOMIAL AT XP,YP
  500 CONTINUE
      XR= XP -X(3)
      YR= YP -Y(3)
      CALL C2HORN (XR,YR,ZP)
C
      IOUTS= 0
      GOTO 2000
C
C     POINT IS OUTSIDE TRIANGULATION
C
 1000 IOUTS= 1
      ZP= 0.
      IF (I3.NE.0) IST= I3
      IF (I2.NE.0) IST= I2
      IF (I1.NE.0) IST= I1
 3000 RETURN
C
C     SAVE VALUES OF I1,I2,I3
 2000 CONTINUE
      I1OLD= I1
      I2OLD= I2
      I3OLD= I3
C
      RETURN
      END
      SUBROUTINE C2SIDE (X,Y,Z,PDL,NSIDES)
C
C     PART OF C2-INTERPOLATION PACKAGE
C
C     COMPUTATION OF COEFFICIENTS FOR POLYNOMIALS ALONG SIDES
C
C     AUTHOR      : A. PREUSSER, 1989
C                   FRITZ-HABER-INSTITUT
C                   DER MAX-PLANCK-GESELLSCHAFT
C                   FARADAYWEG 4-6
C                   D-1000 BERLIN 33
C
C
C     X        X-COORDINATES OF VERTICES
C     Y        Y-COORDINATES OF VERTICES
C     PDL(J,K) PARTIAL DERIVATIVES AT VERTEX K
C                             ZX,ZY
C                         ZXX, ZXY, ZYY
C                    ZXXX, ZXXY, ZXYY, ZYYY
C               ZXXXX, ZXXXY, ZXXYY, ZXYYY, ZYYYY
C              IN THAT ORDER.
C
C     NSIDES   NUMBER OF POLYNOMIALS TO BE DETERMINED.
C              TWO VALUES ARE POSSIBLE:
C              2, FOR SIDE 1 AND 2
C              3, FOR SIDE 1, 2, AND 3
C
C                                  2
C                                  *
C                                /  .
C                              /     .
C                            /        .
C                          /           .
C                        /              .
C             SIDE(1)  /                 .
C                    /                    . SIDE(3)
C                  /                       .
C                /                          .
C              /                             .
C            /                                .
C          /                                   .
C        /                                      .
C    3 *-----------------------------------------* 1
C                      SIDE(2)
C
C
*D    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL X,Y,Z,PDL
C
      COMMON /C2PCO/   Q00(3),Q01(3),Q02(3),Q03(3),Q04(3),Q05(3)
     1,                Q06(3),Q07(3),Q08(3),Q09(3)
     2,         Q12(2),Q13(2),Q14(2),Q15(2),Q16(2),Q17(2),Q18(2)
     3,         Q23(2),Q24(2),Q25(2),Q26(2),Q27(2)
     4,         Q34(2),Q35(2),Q36(2)
     5,         Q45(2),Q11,Q22,Q33,Q44
C
C     Q00...Q09   COEFFICIENTS OF POLYNOMIALS ALONG THE SIDES
C     Q12...Q44   COEFFICIENTS FOR BIVARIATE POLYNOMIAL INSIDE TRIANGLE
C
      COMMON /C2CCO/ DX(3),DY(3),DX2(3),DY2(3),SL2(3)
     1,              AP,BP,CP,DP
C
      COMMON /C2DCO/  Z10(3,3),Z20(3,3),Z30(3,3),Z40(3,3)
     1,               Z21(2,3),Z31(2,3),Z11(3),Z22(3)
C
C
C     OUTPUT OF THIS ROUTINE
C     ----------------------
C     DX,DY     COORDINATE DIFFERENCES OF ENDPOINTS OF SIDES
C     DX2,DY2   DX**2,DY**2
C     SL2       SIDE LENGTH **2
C     ZI0       PARTIAL DERIVATIVES WITH RESPECT TO TRIANGULAR
C               COORDINATES L1,L2,L3 (CORRESPOND TO DIRECTIONS
C               OF SIDES 2,1,3 RESPECTIVELY)
C               FIRST INDEX: VARIABLE
C               SECOND INDEX: POINT
C               I= ORDER OF DERIVATIVE
C     Q00...Q09 COEFFICIENTS FOR POLYNOMIALS ALONG THE THREE
C               SIDES
C
C
      DIMENSION PDL(14,3),X(3),Y(3),Z(3)
C
C     COMPUTE COORDINATE DIFFERENCES FOR SIDES
      DO 60 I=1,3
        NP1= 3
        NP2= 2
        IF (I.EQ.3) NP1= 1
        IF (I.EQ.2) NP2= 1
        XNP1= X(NP1)
        YNP1= Y(NP1)
        DX(I)= X(NP2) - XNP1
        DY(I)= Y(NP2) - YNP1
        DX2(I)= DX(I)*DX(I)
        DY2(I)= DY(I)*DY(I)
        SL2(I)= DX2(I) + DY2(I)
   60 CONTINUE
C
C
C     CONVERTS THE PARTIAL DERIVATIVES AT THE VERTICES OF THE
C     TRIANGLE FOR THE L1-L2 COORDINATE SYSTEM,
C     APPLYING CHAIN RULE
C
C     LOOP ON VARIABLES I
      DO 260 I=1,NSIDES
        J= 3 - MOD(I,3)
        E2XY= 2.*DX(J)*DY(J)
        DX3 = DX2(J)*DX(J)
        E3X2Y= 3.*DX2(J)*DY(J)
        E3Y2X= 3.*DY2(J)*DX(J)
        DY3 = DY2(J)*DY(J)
        DX4 = DX3*DX(J)
        DY4 = DY3*DY(J)
        D4X3Y= 4.*DX3*DY(J)
        D4Y3X= 4.*DY3*DX(J)
        D6X2Y2= 6.*DX2(J)*DY2(J)
C
C       COMPUTE DERIVATIVES IN DIRECTIONS OF SIDES J
C       LOOP ON POINTS K
        NK= 3
        IF (I.GT.2) NK= 2
        DO 230 K=1,NK
          Z10(I,K)= DX(J)*PDL(1,K)  + DY(J)*PDL(2,K)
          Z20(I,K)= DX2(J)*PDL(3,K) + E2XY*PDL(4,K) + DY2(J)*PDL(5,K)
          Z30(I,K)= DX3*PDL(6,K) + E3X2Y*PDL(7,K) +
     .              DY3*PDL(9,K) + E3Y2X*PDL(8,K)
          Z40(I,K)= DX4*PDL(10,K) + D4X3Y*PDL(11,K) +
     .              DY4*PDL(14,K) + D4Y3X*PDL(13,K) +
     .              D6X2Y2*PDL(12,K)
  230   CONTINUE
C
  260 CONTINUE
C
C
C     COMPUTE COEFFICIENTS OF POLYNOMIALS ALONG THE THREE SIDES
C     OF THE TRIANGLE
      DO 300 I=1,NSIDES
        J= 3-MOD(I,3)
C
        NP1= 3
        NP2= 2
        IF (I.EQ.3) NP1=1
        IF (I.EQ.2) NP2=1
C
C
        Q00(I)= Z(NP1)
        Q01(I)= Z10(J,NP1)
        Q02(I)= Z20(J,NP1)*0.5
        Q03(I)= Z30(J,NP1)/6.
        Q04(I)= Z40(J,NP1)/24.
        Z1MZ2= Q00(I) - Z(NP2)
        HQ4= Q04(I)*5.
        HQ42= HQ4 + HQ4
        Q05(I)= Z40(J,NP2)/24. - Z30(J,NP2) +10.5*Z20(J,NP2)
     1          -56.*Z10(J,NP2)  - HQ4
     2          -15.*Q03(I) - 35.*Q02(I) - 70.*Q01(I) -126.*Z1MZ2
        Q06(I)= (-Z40(J,NP2)+23.*Z30(J,NP2))/6.-38.5*Z20(J,NP2)
     1          +196.*Z10(J,NP2) + HQ42
     2          +40.*Q03(I) +105.*Q02(I) +224.*Q01(I) +420.*Z1MZ2
        Q07(I)= 0.25*Z40(J,NP2) -5.5*Z30(J,NP2) +53.*Z20(J,NP2)
     1          -260.*Z10(J,NP2) - HQ42
     2          -45.*Q03(I) -126.*Q02(I) -280.*Q01(I) -540*Z1MZ2
        Q08(I)= -Z40(J,NP2)/6.+3.5*Z30(J,NP2) -32.5*Z20(J,NP2)
     1          +155.*Z10(J,NP2) + HQ4
     2          +24.*Q03(I) +70.*Q02(I) +160.*Q01(I) +315.*Z1MZ2
        Q09(I)= (Z40(J,NP2) -20*Z30(J,NP2))/24. +7.5*Z20(J,NP2)
     1          -35.*Z10(J,NP2) - Q04(I)
     2          -5.*Q03(I) -15.*Q02(I) -35.*Q01(I) -70.*Z1MZ2
C
  300 CONTINUE
C
      RETURN
      END
      SUBROUTINE C2INSD (PDL)
C
C     PART OF C2-INTERPOLATION PACKAGE
C
C     COMPUTATION OF POLYNOMIAL COEFFICIENTS Q12...Q44 FOR
C     BIVARIATE POLYNOMIAL INSIDE TRIANGLE
C     AND
C     CONSTANTS FOR TRANSFORMATION OF
C     CARTESIAN TO TRIANGULAR COORDINATES
C
C     AUTHOR      : A. PREUSSER, 1989
C                   FRITZ-HABER-INSTITUT
C                   DER MAX-PLANCK-GESELLSCHAFT
C                   FARADAYWEG 4-6
C                   D-1000 BERLIN 33
C
*D    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL PDL
      DIMENSION PDL(14,3)
C
C     PDL(J,K) PARTIAL DERIVATIVES AT VERTEX K
C                             ZX,ZY
C                         ZXX, ZXY, ZYY
C                    ZXXX, ZXXY, ZXYY, ZYYY
C               ZXXXX, ZXXXY, ZXXYY, ZXYYY, ZYYYY
C              IN THAT ORDER.
C     NOTATION
C     --------
C               COEFFICIENTS:  QIJ(1)= QIJ
C                              QIJ(2)= QJI
C               PARTIALS    :  ZIJ(1,K)= ZIJ(K)
C                              ZIJ(2,K)= ZJI(K)
C                              ZIJ(K) : Z DIFF I-TIMES TO L1
C                                              J-TIMES TO L2
C                                       AT VERTEX K
C
C
C
      COMMON /C2PCO/  Q00(3),Q01(3),Q02(3),Q03(3),Q04(3),Q05(3)
     1,                Q06(3),Q07(3),Q08(3),Q09(3)
     2,         Q12(2),Q13(2),Q14(2),Q15(2),Q16(2),Q17(2),Q18(2)
     3,         Q23(2),Q24(2),Q25(2),Q26(2),Q27(2)
     4,         Q34(2),Q35(2),Q36(2)
     5,         Q45(2),Q11,Q22,Q33,Q44
C     Q00...Q09   COEFFICIENTS OF POLYNOMIALS ALONG THE SIDES
C     Q12...Q44   COEFFICIENTS FOR BIVARIATE POLYNOMIAL INSIDE TRIANGLE
C
C
C
      COMMON /C2CCO/ DX(3),DY(3),DX2(3),DY2(3),SL2(3)
     1,              AP,BP,CP,DP
C
      COMMON /C2DCO/  Z10(3,3),Z20(3,3),Z30(3,3),Z40(3,3)
     1,               Z21(2,3),Z31(2,3),Z11(3),Z22(3)
C
      PARAMETER (ONE=1.,TWO= 2.,FIVE= 5.)
C
C
C     INPUT TO THIS ROUTINE
C     ---------------------
C     VARIABLES FROM /C2CCO/ AND /C2DCO/ COMPUTED IN C2SIDE
C
C
C     OUTPUT OF THIS ROUTINE
C     ----------------------
C     QIJ  POLYNOMIAL COEFFICIENTS FOR REPRESENTATION INSIDE TRIANGLE
C     AP,BP,CP,DP    CONST. FOR TRANSF. CART. TO TRIANG. SYSTEM
C
C
C--------------------------------------------------------------------
C
      DIMENSION XJ2XI(2),YJ2YI(2)
C
C     COMPUTE CONST. FOR TRANSF. CART.-TRIANG. COORD. SYSTEM
      AD=  DX(2)*DY(1)
      BC=  DX(1)*DY(2)
      DLT= AD-BC
      AP=  DY(1)/DLT
      BP= -DX(1)/DLT
      CP= -DY(2)/DLT
      DP=  DX(2)/DLT
C
C     SOME ADDITIONAL CONSTANTS
      ADBC= AD+BC
      AB= DX(2)*DX(1)
      CD= DY(1)*DY(2)
      ABCD= AB+CD
C
C
C       *****************************
C       *  COMPUTE MIXED PARTIALS   *
C       *****************************
C
      DO 100 I=1,2
        J= 3 - I
        XJ2XI(I)= DX2(J)*DX(I)
        YJ2YI(I)= DY2(J)*DY(I)
        XJ2YI= DX2(J)*DY(I)
        YJ2XI= DY2(J)*DX(I)
        E2XXY= AB*DY(J)*2. + XJ2YI
        E2YYX= CD*DX(J)*2. + YJ2XI
        XJ3XI= XJ2XI(I)*DX(J)
        YJ3YI= YJ2YI(I)*DY(J)
        E3X2XY= 3.*XJ2XI(I)*DY(J) + XJ2YI*DX(J)
        E3Y2YX= 3.*YJ2YI(I)*DX(J) + YJ2XI*DY(J)
        E3J2= 3.*(XJ2YI*DY(J) + YJ2XI*DX(J))
        DO 40 K=1,3
          Z21(I,K)= XJ2XI(I)*PDL(6,K) + YJ2YI(I)*PDL(9,K) +
     .              E2XXY*PDL(7,K) + E2YYX*PDL(8,K)
          Z31(I,K)= XJ3XI*PDL(10,K) + YJ3YI*PDL(14,K) +
     .              E3X2XY*PDL(11,K) + E3Y2YX*PDL(13,K) +
     .              E3J2*PDL(12,K)
   40   CONTINUE
  100 CONTINUE
C
      DX1X22= DX2(1)*DX2(2)
      DY1Y22= DY2(1)*DY2(2)
      E2IJX= 2.*(XJ2XI(1)*DY(1) + XJ2XI(2)*DY(2))
      E2IJY= 2.*(YJ2YI(1)*DX(1) + YJ2YI(2)*DX(2))
      E12  = DX2(2)*DY2(1) + DX2(1)*DY2(2) + 4.*AB*CD
      DO 200 K=1,3
        Z11(K)= AB*PDL(3,K) + ADBC*PDL(4,K) + CD*PDL(5,K)
        Z22(K)= DX1X22*PDL(10,K) + DY1Y22*PDL(14,K) +
     .          E2IJX*PDL(11,K) + E2IJY*PDL(13,K) +
     .          E12  *PDL(12,K)
  200 CONTINUE
C
C
C        ***************************************
C        *  COMPUTE Q1I- AND Q2I-COEFFICIENTS  *
C        ***************************************
C
      Q11= Z11(3)
      Q22= Z22(3)*0.25
      T6= ONE/6.
      DO 500 I=1,2
        J= 3-I
C
        B= SL2(I)
        D= ABCD
        BDNO= (ABS(B)+ABS(D))*0.5
        B= B/BDNO
        D= D/BDNO
        DB= D/B
        DB2= DB*DB
C
        Q18(I)= 9.*Q09(I)*DB
        Q12(I)= Z21(J,3)*0.5
        Q13(I)= Z31(J,3)*T6
C
        Q184= Q18(I)*4.
        Q134= Q13(I)*4.
        Z31JJ6= Z31(J,J)*T6
        Z31JJ2= Z31(J,J)*0.5
        ZTDIF= Z10(I,3)-Z10(I,J)
        Q14(I)= -Z31JJ6 +2.5*Z21(J,J) -15.*Z11(J)  +Q18(I)
     1         -Q134      -10.*Q12(I) -20.*Q11 -35.*ZTDIF
*
        Q15(I)= Z31JJ2   -7.*Z21(J,J) +39.*Z11(J) - Q184
     1         +6.*Q13(I) +20.*Q12(I) +45.*Q11 +84.*ZTDIF
*
        Q16(I)= -Z31JJ2 +6.5*Z21(J,J) -34.*Z11(J)+ 6.*Q18(I)
     1         -Q134      -15.*Q12(I) -36.*Q11 -70.*ZTDIF
*
        Q17(I)= Z31JJ6   -2.*Z21(J,J) +10.*Z11(J) - Q184
     1         +   Q13(I) + 4.*Q12(I) +10.*Q11 +20.*ZTDIF
C
C
        Q26(I)=  7.*Q17(I)*DB - 28.*Q08(I)*DB2
        Q27(I)=  8.*Q18(I)*DB - 36.*Q09(I)*DB2
C
        Z20JM3= Z20(I,J) - Z20(I,3)
        Z224= Z22(J)*0.25
        Q263= Q26(I)*3.
        Q223= Q22*3.
        Q25(I)= -6.*Q27(I) -Q263 - Q22
     1          -1.5*Z21(I,3) +Z224
     2          -1.5*Z21(I,J) +3.*Z20JM3
C
        Q24(I)=  8.*Q27(I) +Q263 + Q223
     1          +4.*Z21(I,3) - 0.5*Z22(J)
     2          +3.5*Z21(I,J) -7.5*Z20JM3
C
        Q23(I)= -3.*Q27(I) -Q26(I) -Q223
     1          -3.*Z21(I,3) + Z224
     2          -2.*Z21(I,J) +5.*Z20JM3
 500  CONTINUE
C
C
C       *****************************************
C       *   COMPUTE Q3I- AND Q4I-COEFFICIENTS   *
C       *             (LAST TEN)                *
C       *****************************************
C
      B=  DX(1)*DX(3) + DY(1)*DY(3)
      D= -DX(2)*DX(3) - DY(2)*DY(3)
      BDNO= (ABS(B)+ABS(D))*0.5
      B= B/BDNO
      D= D/BDNO
      BB= B*B
      DD= D*D
      BBPDD= BB + DD
      BPD= B + D
      BPD2= BPD*BPD
      BMD= B - D
      BMD2= BMD*BMD
      BD= B*D
C
      T3= ONE/3.
      T23= TWO/3.
C
      BB9= BB*9.
      BB11= BB*11.
C
      DD2= DD*2.
      DD18= DD*18.
C
      BD3= BD*3.
      BD9= BD*9.
      BD15= BD*15.
C
      C1=  DD*32. - BD15 + BB*25.
      C2=  DD*10. + BD15 + BB*17.
      C3=  DD*12. - BD3  + BB*5.
      C4=  DD2    + BD3  + BB9
      C5=  DD*25. - BD9  + BB*18.
      C6=  DD*6.  + BD15 + BB*13.
      C7=  DD18   - BD3  + BB11
      C8=  DD2    + BD15 + BB9
      C9=  DD18   + BD9  + BB11
      C0=  DD*8.  + BD15 + BB*15.
C
      PART1= (Q03(1)*C1 + Q03(2)*C2 +
     .       Q04(1)*C3 + Q04(2)*C4 +
     .       Q13(1)*C5 + Q13(2)*C6 +
     .       Q14(1)*C3 + Q14(2)*C4 +
     .       Q23(1)*C7 + Q23(2)*C8 +
     .       Q24(1)*C9 + Q24(2)*C0)*3.
      PART2= Z31(2,1)*(3.5*BBPDD  - BD3) +
     .       Z31(1,2)*BBPDD*2. +
     .       Z40(1,2)*(-DD*0.25 -BD*0.375 - BB*1.125) +
     .       Z40(2,1)*(-DD*1.5  +BD*0.375 - BB*0.625) -
     .       (Z30(1,2)*C2 + Z30(2,1)*C1)*0.5
      PART3= ( 14.*Q27(2)  +  27.*Q26(2) + 30.*Q25(2)
     .         -T3*Q18(2)  -  84.*Q27(1) - 15.*Q26(1)
     .        +30.*Q25(1)  + (226.+T23)*Q18(1) + 63.*Q17(1)
     .       -399.*Q09(1)  -  84.*Q08(1))*DD
      PART4=  (151.*Q27(2)    +  96.*Q26(2) + 60.*Q25(2)
     .       -(81.+T23)*Q18(2) -  21.*Q17(2)  - 3.*Q09(2)
     .       + 129.*Q27(1)    +  96.*Q26(1) + 60.*Q25(1)
     .       -(53.+T23)*Q18(1) -  21.*Q17(1) - 39.*Q09(1))*BD
      PART5= (-  91.*Q27(2)    - 15.*Q26(2) +  30.*Q25(2)
     .       +(230.+T23)*Q18(2) + 63.*Q17(2) - 399.*Q09(2)
     .       -      84.*Q08(2) + 21.*Q27(1) + 27.*Q26(1)
     .       +      30.*Q25(1) - (4.+T3)*Q18(1))*BB
      Q45(1)= -(PART1+PART2+PART3+PART4+PART5)/BPD2
C
      Q45(2)=  Q45(1) +
     .        ((-36.*(Q09(1)-Q09(2)) + 20.*(Q18(1)-Q18(2))
     .          - 8.*(Q27(1)-Q27(2)))*BD
     .      +(45.*(Q03(1)-Q03(2)) +  9.*(Q04(1)-Q04(2))
     .      + 36.*(Q13(1)-Q13(2)) +  9.*(Q14(1)-Q14(2))
     .      + 27.*(Q23(1)-Q23(2)) +  9.*(Q24(1)-Q24(2))
     .      +  7.*(Q27(1)-Q27(2)) -  4.*(Q18(1)-Q18(2))
     .      + 0.375*(Z40(1,2)-Z40(2,1)) + 7.5*(Z30(1,2)-Z30(2,1))
     .      - 1.5*(Z31(1,2) - Z31(2,1)))* BMD2)/BPD2
      DM2B= D-B*2
C
      PART1= ((Q03(2)-Q03(1))*5. + (Q13(2)-Q13(1))*4. +
     .        (Q23(2)-Q23(1))*3. +  Q24(2)-Q24(1)     +
     .         Q14(2)-Q14(1)     +  Q04(2)-Q04(1))*DM2B*24.
      PART2= ((Z30(2,1)-Z30(1,2))*20. + (Z31(1,2)-Z31(2,1))*4. +
     .         Z40(2,1)-Z40(1,2))*DM2B
      PART3= (2.*Q27(2)-Q18(2)-7.*Q27(1)+8.*Q18(1)-9.*Q09(1))*D*8.
      PART4= (2.*Q27(1)-Q18(1)-7.*Q27(2)+8.*Q18(2)-9.*Q09(2))*B*8.
      PART5= (Q45(2)*8.-Q45(1)*16.)*BPD
      Q36(1)= -(PART1+PART2+PART3+PART4+PART5)/(24.*BPD)
C
      T24= ONE/24.
      T56= FIVE/6.
      Q36(2)=  - Z40(2,1)*T24 + Z31(2,1)*T6
     .         - Z30(2,1)*T56 + Q45(2) - Q24(2)
     .         - 3.*Q23(2) - Q14(2) - 4.*Q13(2) - Q04(2)
     .         - 5.*Q03(2) + Z40(1,2)*T24 - Z31(1,2)*T6
     .         + Z30(1, 2)*T56 - Q45(1) + Q36(1)
     .         + Q24(1) + 3.*Q23(1) + Q14(1) + 4.*Q13(1)
     .         + Q04(1) +  5.*Q03(1)
C
      PART1= (5.*Z40(2,1)-14.*Z31(2,1)+7.*Z40(1,2))/36.
      PART2= 10.*(-Q26(2)-Q25(2)-Q26(1)-Q25(1)) -
     .        8.*Q23(2)    -T23*Q36(1) -  6.*Q23(1)
      PART3= (8.*Z30(2,1)  + 5.*Q45(2) - 35.*Q27(2)
     .        -48.*(Q24(2) + Q13(2) + Q03(1))
     .        -21.*(Q14(2) + Q04(2)) -60.*Q03(2)
     .        -2.*Z31(1,2) + 10.*Z30(1,2)
     .        -Q45(1)      - 35.*Q27(1)  -42.*Q24(1)
     .        -15.*Q14(1)  - 37.5*Q13(1) -15.*Q04(1))/4.5
      Q35(2)= -PART1-PART2-PART3
C
      PART1= 2.*(Q36(2)-Q36(1)+Q23(2)-Q23(1))
      PART2= 3.*(Q13(2)-Q13(1))
      PART3= 4.*(Q03(2)-Q03(1))
      PART4= T23*(Z30(2,1)-Z30(1,2))
      PART5= T6* (Z31(1,2)-Z31(2,1))
      Q35(1)= Q35(2) +
     .        PART1 + PART2 + PART3 + PART4 + PART5
C
      Q33= -  2.*Q23(2) -  3.*Q13(2) - 4.*Q03(2)
     .     - Z31(1,2)*T6 + T23*Z30(1,2) + 2.*Q36(1) + Q35(1)
C
      Q34(2)=  -0.75*Q33 + Z31(2,1)*T24
     .         - 1.5*Q36(2) -1.25*Q35(2) -0.5*Q23(1) -0.25*Q13(1)
C
      Q34(1)=  -0.75*Q33 -0.5*Q23(2)
     .         -0.25*Q13(2)+Z31(1,2)*T24 -1.5*Q36(1) -1.25*Q35(1)
C
      Q44= + Z40(2,1)*T24 - Q45(2)
     .     - Q23(2) - 2.*Q13(2) - 3.*Q03(2) -   Z31(1, 2)*T6
     .     + 0.5*Z30(1,2) + 3.*Q36(1) + 2.*Q35(1)
     .     - Q24(1) -  Q14(1) -  Q04(1)
C
      RETURN
      END
      SUBROUTINE C2HORN (X,Y,Z)
C
C     PART OF C2-INTERPOLATION PACKAGE
C
C     EVALUATION OF BIVARIATE NONIC POLYNOMIAL
C
C     AUTHOR      : A. PREUSSER, 1989
C                   FRITZ-HABER-INSTITUT
C                   DER MAX-PLANCK-GESELLSCHAFT
C                   FARADAYWEG 4-6
C                   D-1000 BERLIN 33
C
C
C    INPUT PARAMETERS
C    ----------------
C    X,Y      CARTESIAN COORDINATES RELATIVE TO VERTEX (3)
C             OF TRIANGLE
C
C    OUTPUT PARAMETER
C    ----------------
C    Z        VALUE OF POLYNOMIAL AT X,Y
C
*D    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL X,Y
C
      COMMON /C2PCO/  Q00(3),Q01(3),Q02(3),Q03(3),Q04(3),Q05(3)
     1,                Q06(3),Q07(3),Q08(3),Q09(3)
     2,         Q12(2),Q13(2),Q14(2),Q15(2),Q16(2),Q17(2),Q18(2)
     3,         Q23(2),Q24(2),Q25(2),Q26(2),Q27(2)
     4,         Q34(2),Q35(2),Q36(2)
     5,         Q45(2),Q11,Q22,Q33,Q44
C
      COMMON /C2CCO/ DX(3),DY(3),DX2(3),DY2(3),SL2(3)
     1,              AP,BP,CP,DP
C
C     INPUT TO THIS ROUTINE
C     ---------------------
C     ALL VARIABLES FROM /C2PCO/ (POLYNOMIAL COEFFICIENTS)
C     AP,BP,CP,DP    CONSTANTS DERIVED FROM COORDINATE DIFFERENCES
C                    ON SIDES
C
C ------------------------------------------------------------------
C
C     TRANSFORMATION FROM CARTESIAN TO TRIANGULAR COORDINATES
C                                                (U=L1, V=L2)
      U= AP*X + BP*Y
      V= CP*X + DP*Y
C
C     EVALUATION BY HORNER'S SCHEME (BIVARIATE)
      H0= Q00(2)+V*(Q01(1)+V*(Q02(1)+V*(Q03(1)+V*(Q04(1)+V*(Q05(1)
     1          +V*(Q06(1)+V*(Q07(1)+V*(Q08(1)+V*Q09(1)))))))))
      H1= Q01(2) +V*(Q11    +V*(Q12(1) +V*(Q13(1) +V*(Q14(1)
     1           +V*(Q15(1) +V*(Q16(1) +V*(Q17(1) +V*Q18(1))))))))
      H2= Q02(2) +V*(Q12(2) +V*(Q22    +V*(Q23(1) +V*(Q24(1)
     1           +V*(Q25(1) +V*(Q26(1) +V*Q27(1)))))))
      H3= Q03(2) +V*(Q13(2) +V*(Q23(2) +V*(Q33    +V*(Q34(1)
     1           +V*(Q35(1) +V*Q36(1))))))
      H4= Q04(2) +V*(Q14(2) +V*(Q24(2) +V*(Q34(2) +V*(Q44
     1           +V*Q45(1)))))
      H5= Q05(2) +V*(Q15(2) +V*(Q25(2) +V*(Q35(2) +V*Q45(2))))
      H6= Q06(2) +V*(Q16(2) +V*(Q26(2) +V*Q36(2)))
      H7= Q07(2) +V*(Q17(2) +V*Q27(2))
      H8= Q08(2) +V*Q18(2)
      Z = H0 +U*(H1 +U*(H2 +U*(H3 +U*(H4
     1       +U*(H5 +U*(H6 +U*(H7 +U*(H8 +U*Q09(2)))))))))
C
      RETURN
      END
      SUBROUTINE C1GRID (ND,XD,YD,ZD,NX,NY,XX,YY,Z,NDIMX,WK,NW,
     1                   IWK,IW)
C
      DIMENSION XD(ND),YD(ND),ZD(ND)
      DIMENSION XX(NX),YY(NY)
      DIMENSION Z(NDIMX,NY)
      DIMENSION WK(ND,6),IWK(IW)
C
C**********************************************************************
C
C     COMPUTATION OF A RECTANGULAR GRID FROM SCATTERED DATA
C     ===========================================================
C     USING A C1 QUINTIC POLYNOMIAL REPRESENTATION OVER TRIANGLES
C
C
C     AUTHOR      : A. PREUSSER, 1989
C                   FRITZ-HABER-INSTITUT
C                   DER MAX-PLANCK-GESELLSCHAFT
C                   FARADAYWEG 4-6
C                   D-1000 BERLIN 33
C
C     INPUT PARAMETERS
C     ----------------
C
C     ND          NUMBER OF DATA POINTS
C     XD,YD,ZD    ARRAYS OF LENGTH ND  WITH
C                 X,Y,Z-COORDINATES OF DATA POINTS
C     NX,NY       NUMBER OF GRID LINES IN X- AND
C                 Y-DIRECTION.
C     XX,YY       ARRAYS OF LENGTH NX, RESP. NY,
C                 WITH COORDINATES FOR GRID LINES
C                 X= XX(I), Y= YY(J), I=1,NX; J=1,NY,
C                 IN ASCENDING ORDER.
C     Z           OUTPUT PARAMETER, SEE BELOW.
C     NDIMX       FIRST DIMENSION OF THE TWO-DIMENSIONAL
C                 ARRAY Z AS DECLARED IN THE CALLING
C                 (SUB-)PROGRAM.
C     WK          REAL WORK ARRAY OF LENGTH NW.
C     NW          DIMENSION OF ARRAY WK,
C                 NW.GE.6*ND.
C     IWK         INTEGER WORK ARRAY OF LENGTH IW.
C     IW          DIMENSION OF ARRAY IWK,
C                 IW.GE.7*ND.
C
C
C     OUTPUT PARAMETER
C     ----------------
C
C     Z           TWO-DIMENSIONAL ARRAY WITH DIMENSION
C                 (NDIMX,NY) CONTAINING THE INTERPOLATED
C                 Z-VALUES AT THE INTERSECTIONS OF THE
C                 GRID LINES (=GRID POINTS).
C                 Z(I,J) CORRESPONDS TO THE POINT WITH
C                 THE COORDINATES XX(I),YY(J).
C                 IF XX(I),YY(J) IS OUTSIDE THE TRIANGULATION
C                 OF THE ND DATA POINTS, Z(I,J) IS NOT ALTERED.
C
C     METHODS
C     -------
C
C     A CONVEX TRIANGULAR MESH IS FORMED IN THE X-Y PLANE,
C     WITH THE PROJECTIONS OF THE ND DATA POINTS AS NODES.
C     (SUBROUTINE TRMESH FROM ACM ALG. 624 BY R.RENKA)
C
C     PARTIAL DERIVATIVES ZX,ZY,ZXX,ZXY,ZYY ARE ESTIMATED
C     AT THE NODES.
C     (SUBSEQUENT CALLS TO GRANKA, A MODIFIED VERSION OF GRADG
C      FROM ACM ALG. 624 BY R.RENKA)
C
C     A BIVARIATE QUINTIC HERMITE POLYNOMIAL IS FORMED FOR
C     EACH TRIANGLE.
C     (ACM ALG. 526 AND 626 BY H.AKIMA AND A.PREUSSER)
C
C     THE POLYNOMIAL IS EVALUATED BY HORNER'S METHOD, IN CASE
C     THE REQUESTED POINT XX(I),YY(J) IS INSIDE THE TRIANGULATION.
C     IF THE POINT IS OUTSIDE, Z(I,J) IS NOT ALTERED.
C**********************************************************************
C
      DIMENSION XL(3),YL(3),ZL(3),PDL(5,3)
C
C
C     ERROR CHECKS
C
      IF (ND.LT.3) THEN
        WRITE (*,901) ND
        GOTO 1000
       END IF
      IF (NX.LT.1) THEN
        WRITE (*,902) NX
        GOTO 1000
       END IF
      IF (NY.LT.1) THEN
        WRITE (*,903) NY
        GOTO 1000
       END IF
      IF (NX.GT.NDIMX) THEN
        WRITE (*,904) NX,NDIMX
        GOTO 1000
       END IF
      IF (NW.LT.6*ND)  THEN
        WRITE (*,905) NW,6*ND
        GOTO 1000
       END IF
      IF (IW.LT.7*ND)  THEN
        WRITE (*,906) IW,7*ND
        GOTO 1000
       END IF
      IXERR= 0
      DO 20 I=2,NX
        IF (XX(I)-XX(I-1).GT.0.) GOTO 20
        IXERR= 1
   20 CONTINUE
      IF (IXERR.EQ.1)  WRITE (*,907) I,XX(I),I-1,XX(I-1)
      IYERR= 0
      DO 30 I=2,NY
        IF (YY(I)-YY(I-1).GT.0.) GOTO 30
        IYERR= 1
   30 CONTINUE
      IF (IYERR.EQ.1)  WRITE (*,908) I,YY(I),I-1,YY(I-1)
  901 FORMAT (' ***ERROR*** ND.LT.3, ND=',I10)
  902 FORMAT (' ***ERROR*** NX.LT.1, NX=',I10)
  903 FORMAT (' ***ERROR*** NY.LT.1, NY=',I10)
  904 FORMAT (' ***ERROR*** NX.GT.NXDIM. NX,NXDIM=',I10,',',I10)
  905 FORMAT (' ***ERROR*** NW.LT.6*ND. NW,6*ND=',I10,',',I10)
  906 FORMAT (' ***ERROR*** IW.LT.7*ND. IW,7*ND=',I10,',',I10)
  907 FORMAT ('0***ERROR***          '/
     A        '      X-COORDINATE', I5,'=',E15.7/
     B        ' .LE. X-COORDINATE', I5,'=',E15.7)
  908 FORMAT ('0***ERROR***          '/
     A        '      Y-COORDINATE', I5,'=',E15.7/
     B        ' .LE. Y-COORDINATE', I5,'=',E15.7)
C
C     TRIANGULATION AND COMPUTATION OF PARTIAL DERIVATIVES
      N2= 6*ND
      CALL C1INI (ND,XD,YD,ZD,WK,NW,IWK,IW)
C
C
C COMPUTE INTERPOLATED VALUES  Z(I,J)
C
C     LOOPS TRIANGLES
      INDF= 1
      N21= N2-1
      DO 500 K1=1,ND-2
        I1= K1
        INDL= IWK(N21+K1)
        DO 400 INDX=INDF,INDL
          I2= IWK(INDX)
          I3= IWK(INDX+1)
          IF (INDX.EQ.INDL) I3= IWK(INDF)
          IF (I2.LT.K1 .OR. I3.LT.K1) GOTO 400
C
C         START COMPUTATIONS FOR TRIANGLE
C
C         EXTREME VALUES OF VERTEX COORDINATES
          XMIN= AMIN1(XD(I1),XD(I2),XD(I3))
          XMAX= AMAX1(XD(I1),XD(I2),XD(I3))
          YMIN= AMIN1(YD(I1),YD(I2),YD(I3))
          YMAX= AMAX1(YD(I1),YD(I2),YD(I3))
C
C         FIND MIN AND MAX INDICES IN X
C
          IXMIN= 0
          DO 110 I=1,NX
            IF (XX(I).LT.XMIN) GOTO 110
            IXMIN= I
            GOTO 115
  110     CONTINUE
          IF (IXMIN.EQ.0) GOTO 400
  115     CONTINUE
C
          IXMAX= 0
          DO 120 I=NX,1,-1
            IF (XX(I).GT.XMAX) GOTO 120
            IXMAX= I
            GOTO 125
  120     CONTINUE
          IF (IXMAX.EQ.0) GOTO 400
  125     CONTINUE
C
C
C         LOOPS GRID POINTS
          ITFLAG= 0
          DO 390 J=1,NY
            IF (YY(J).LT.YMIN .OR. YY(J).GT.YMAX) GOTO 390
            DO 370 I=IXMIN,IXMAX
C
C             START COMPUTATIONS FOR GRID POINT
C
C             CHECK IF GRID POINT INSIDE TRIANGLE
              IF ((XD(I1)-XX(I))*(YD(I2)-YY(J)) -
     1            (YD(I1)-YY(J))*(XD(I2)-XX(I)) .LT. 0. .OR.
     2            (XD(I2)-XX(I))*(YD(I3)-YY(J)) -
     3            (YD(I2)-YY(J))*(XD(I3)-XX(I)) .LT. 0. .OR.
     4            (XD(I3)-XX(I))*(YD(I1)-YY(J)) -
     5            (YD(I3)-YY(J))*(XD(I1)-XX(I)) .LT.0) GOTO 370
C
C             CHECK IF COEFFICIENTS ALREADY COMPUTED
              IF (ITFLAG.NE.0) GOTO 350
              ITFLAG= 1
C
C             LOAD COORDINATES
              XL(1)= XD(I1)
              YL(1)= YD(I1)
              XL(2)= XD(I2)
              YL(2)= YD(I2)
              XL(3)= XD(I3)
              YL(3)= YD(I3)
C
              ZL(1)= ZD(I1)
              ZL(2)= ZD(I2)
              ZL(3)= ZD(I3)
C             LOAD PARTIAL DERIVATIVES
              DO 300 K=1,5
                PDL(K,1)= WK(I1,K)
                PDL(K,2)= WK(I2,K)
                PDL(K,3)= WK(I3,K)
  300         CONTINUE
C
C             COMPUTE COEFFICIENTS ALONG SIDES
              CALL C1SIDE (XL,YL,ZL,PDL,2)
C
C             COMPUTE COEFFICIENTS FOR INSIDE OF TRIANGLE
C             AND CONSTANTS FOR TRANSF. CARTESIAN TO TRIANGULAR
              CALL C1INSD (PDL)
C
C             EVALUATE POLYNOMIAL AT XP,YP
  350         CONTINUE
              XR= XX(I) -XL(3)
              YR= YY(J) -YL(3)
              CALL C1HORN (XR,YR,Z(I,J))
C
  370       CONTINUE
  390     CONTINUE
  400   CONTINUE
        INDF= INDL + 1
  500 CONTINUE
C
 1000 RETURN
      END
      SUBROUTINE C1INI (ND,XD,YD,ZD,WK,NW,IREN,IR)
C
      DIMENSION XD(ND),YD(ND),ZD(ND)
      DIMENSION WK(ND,6),IREN(IR)
C
C**********************************************************************
C
C     TRIANGULATE DATA AND COMPUTE PARTIAL DERIVATIVES UP TO 2ND ORDER
C     ================================================================
C
C
C     AUTHOR      : A. PREUSSER, 1989
C                   FRITZ-HABER-INSTITUT
C                   DER MAX-PLANCK-GESELLSCHAFT
C                   FARADAYWEG 4-6
C                   D-1000 BERLIN 33
C
C     INPUT PARAMETERS
C     ----------------
C
C     ND          NUMBER OF DATA POINTS
C     XD,YD,ZD    ARRAYS OF LENGTH ND  WITH
C                 X,Y,Z-COORDINATES OF DATA POINTS
C
C     OUTPUT PARAMETERS
C     -----------------
C
C     WK          REAL WORK ARRAY OF LENGTH NW.
C     NW          DIMENSION OF ARRAY WK,
C                 NW.GE.6*ND.
C                 USAGE OF COLUMNS OF ARRAY WK:
C
C                                    1    2
C                                   ZX    ZY
C
C                                3     4     5
C                               ZXX   ZXY   ZYY
C
C                 THE 6TH COLUMN IS WORK SPACE.
C
C     IREN        INTEGER WORK ARRAY OF LENGTH IR.
C     IR          DIMENSION OF ARRAY IREN,
C                 IR.GE.7*ND.
C
C     METHODS
C     -------
C
C     A CONVEX TRIANGULAR MESH IS FORMED IN THE X-Y PLANE,
C     WITH THE PROJECTIONS OF THE ND DATA POINTS AS NODES.
C     (SUBROUTINE TRMESH FROM ACM ALG. 624 BY R.RENKA)
C
C     PARTIAL DERIVATIVES ZX,ZY,ZXX,ZXY,ZYY ARE ESTIMATED
C     AT THE NODES.
C     (SUBSEQUENT CALLS TO GRANKA, A MODIFIED VERSION OF GRADG
C      FROM ACM ALG. 624 BY R.RENKA)
C
C**********************************************************************
C
C
C
C     ERROR CHECKS
C
      IF (ND.LT.3) THEN
        WRITE (*,901) ND
        GOTO 1000
       END IF
      IF (NW.LT.6*ND)  THEN
        WRITE (*,905) NW,6*ND
        GOTO 1000
       END IF
      IF (IR.LT.7*ND)  THEN
        WRITE (*,906) IR,7*ND
        GOTO 1000
       END IF
  901 FORMAT (' ***ERROR*** ND.LT.3, ND=',I10)
  905 FORMAT (' ***ERROR*** NW.LT.6*ND. NW,6*ND=',I10,',',I10)
  906 FORMAT (' ***ERROR*** IR.LT.7*ND. IR,7*ND=',I10,',',I10)
C
C     TRIANGULATION
      N2= 6*ND
      CALL TRMESH (ND,XD,YD,IREN,IREN(N2),IER)
C     TRMESH IS PART OF ACM ALGORITHM 624 BY R.RENKA
      IF (IER.EQ.3) WRITE (*,911)
      IF (IER.EQ.3) GOTO 1000
  911 FORMAT (' ***ERROR*** ALL POINTS COLLINEAR')
C
C     LOOPS TRIANGLES FOR DUPLICATE POINTS
      INDF= 1
      N21= N2-1
      DO 100 K1=1,ND-1
        IND1= K1
        INDL= IREN(N21+K1)
        DO 40 INDX=INDF,INDL
          IND2= IREN(INDX)
          IF (IND2.LT.K1) GOTO 40
          IF (ABS(XD(IND1)-XD(IND2))+
     .        ABS(YD(IND1)-YD(IND2)).NE.0.) GOTO 40
          WRITE (*,912) IND1,IND2
          GOTO 1000
   40   CONTINUE
        INDF= INDL + 1
  100 CONTINUE
  912 FORMAT (' ***ERROR*** DUPLICATE POINTS',2I10)
C
C
C     COMPUTE PARTIAL DERIVATIVES BY RENKA'S METHOD
C
C     GRANKA IS A MODIFIED VERSION OF GRADG FROM
C     ACM ALGORITHM 624
C
      DO 60 J=1,6
        DO 60 I=1,ND
   60     WK(I,J)= 0.
C
C     COMPUTE ZX AND ZY
C     L= NUMBER OF ITERATIONS
      L= 15
      E= 0.
      CALL GRANKA (ND,XD,YD,ZD,IREN,IREN(N2),E,L,WK,IER)
C
      L= 15
C     COMPUTE ZXX AN ZXY
      CALL GRANKA (ND,XD,YD,WK,IREN,IREN(N2),E,L,WK(1,3),IER)
C
C     COMPUTE ZYX AND ZYY
      L= 15
      CALL GRANKA (ND,XD,YD,WK(1,2),IREN,IREN(N2),E,L,WK(1,5),IER)
C
C
C     COMPUTE/STORE ZXY AND ZYY
C
      DO 80 I=1,ND
        WK(I,4)= (WK(I,4)+WK(I,5))*0.5
        WK(I,5)= WK(I,6)
   80 CONTINUE
C
 1000 CONTINUE
      RETURN
      END
      SUBROUTINE C1PNT (ND,XD,YD,ZD,XP,YP,PD,
     .                  IREN,IST,ZP,IOUTS)
      DIMENSION IREN(7*ND)
      DIMENSION XD(ND), YD(ND), ZD(ND), PD(ND,5)
     A,         X(3),Y(3),Z(3),PDL(5,3)
C
C***********************************************************
C
C     INTERPOLATE Z-VALUE OF A POINT IN A C1-REPRESENTATION
C     -----------------------------------------------------
C                   ON A SET OF TRIANGLES
C                   ---------------------
C
C
C     AUTHOR:
C     A.PREUSSER
C     FRITZ-HABER-INSTITUT DER MPG
C     FARADAYWEG 4-6
C     D-1000 BERLIN 33
C
C
C
C
C GIVEN A TRIANGULATION OF A SET OF POINTS IN THE PLANE,
C IN A DATA STRUCTURE AS DEVELOPED BY R. RENKA FOR ACM ALG. 624,
C AND VALUES AND PARTIAL DERIVATIVES AT THESE POINTS,
C THIS ROUTINE COMPUTES A VALUE ZP AT XP,YP.
C
C THE PREREQUISITE VALUES CAN BE COMPUTED BY SUBROUTINE C1INI.
C
C INPUT PARAMETERS
C ----------------
C       ND -    NUMBER OF NODES IN THE MESH.
C               ND .GE. 3.
C               THE VALUE OF ND MUST BE EXACTLY
C               THE SAME AS SPECIFIED FOR THE CALL
C               OF SUBROUTINE C2INI.
C
C       XD,YD - VECTORS OF COORDINATES OF THE
C               NODES IN THE MESH.
C
C          ZD - VECTOR OF DATA VALUES AT THE
C               NODES.
C
C       XP,YP - COORDINATES OF A POINT AT WHICH
C               FUNCTION VALUE ZP IS REQUIRED
C
C          PD - ND BY 5 ARRAY CONTAINING
C               PARTIAL DERIVATIVES AT THE
C               NODES.
C               FIRST INDEX INDICATES POINT,
C               SECOND INDEX DERIVATIVE:
C               COLUMN 1: ZX
C                      2: ZY
C                      3: ZXX
C                      4: ZXY
C                      5: ZYY
C
C        IREN - INTEGER ARRAY OF LENGTH 7*ND.
C               CAN BE COMPUTED BY A CALL TO SUBROUTINE
C               C2INI. IT CONTAINS THE DATA FOR
C               THE DEFINITION OF THE TRIANGLES.
C               (ARRAYS IADJ AND IEND OF ALG. 624)
C
C         IST - INDEX OF THE STARTING NODE IN
C               THE SEARCH FOR A TRIANGLE CON-
C               TAINING (XP,YP).  1 .LE. IST
C               .LE. ND.  THE OUTPUT VALUE OF
C               IST FROM A PREVIOUS CALL MAY
C               BE A GOOD CHOICE.
C
C     IADJ AND IEND MAY BE CREATED BY TRMESH (ACM ALG 624)
C
C   INPUT PARAMETERS OTHER THAN IST ARE NOT ALTERED BY THIS
C   ROUTINE.
C
C OUTPUT PARAMETERS
C -----------------
C
C         IST - INDEX OF ONE OF THE VERTICES OF
C               THE TRIANGLE CONTAINING (XP,YP),
C
C          ZP - VALUE TO BE COMPUTED AT (XP,YP)
C               =0, IF XP,YP IS OUTSIDE THE TRIANGULATION
C
C       IOUTS - =0, IF (XP,YP) INSIDE THE TRIANGULATION
C                   AND ZP CONTAINS THE INTERPOLATED VALUE.
C               =1, IF (XP,YP) IS OUTSIDE THE TRIANGULATION
C                   AND ZP=0.
C
C     REMARK
C     ------
C     THIS ROUTINE SAVES THE INDICES OF THE THREE VERTICES OF
C     THE TRIANGLE IN WHICH THE POINT XP,YP IS FOUND.
C     IF, IN THE NEXT CALL, THE SAME INDICES ARE FOUND,
C     IT ASSUMES THAT /C1PCO/ AND /C1CCO/ HAVE NOT BEEN
C     ALTERED BY THE USER. THEN THE COEFFICIENTS FOR
C     THE TRIANGLE ARE NOT CALCULATED BUT TAKEN FROM /C1PCO/
C     AND /C1CCO/.
C***********************************************************
C
      SAVE I1OLD,I2OLD,I3OLD,X,Y
      DATA I1OLD,I2OLD,I3OLD /0,0,0/
C
C     ERROR CHECKS
      IF (ND.LT.3) THEN
         WRITE (*,901) ND
         GOTO 3000
       END IF
      IF (IST.LT.1) THEN
         WRITE (*,902) IST
         GOTO 3000
       END IF
      IF (IST.GT.ND) THEN
         WRITE (*,903) IST,ND
         GOTO 3000
       END IF
  901 FORMAT (' ***ERROR*** ND.LT.3, ND= ',I10)
  902 FORMAT (' ***ERROR*** IST.LT.1, IST=',I10)
  903 FORMAT (' ***ERROR*** IST.GT.ND, IST,ND=',I10,',',I10)
C
C     FIND A TRIANGLE CONTAINING P
C
      N2= 6*ND
      CALL TRFIND(IST,XP,YP,XD,YD,IREN,IREN(N2),I1,I2,I3)
      IF (I1 .EQ. 0 .OR. I2.EQ.0 .OR. I3.EQ.0) GOTO 1000
      IST = I1
C
C     CHECK IF SAME TRIANGLE AS BEFORE
C     (COUNTER-CLOCKWISE ORDER OF VERTICES)
      IF (I1.EQ.I1OLD .AND. I2.EQ.I2OLD .AND. I3.EQ.I3OLD
     1.OR.I1.EQ.I2OLD .AND. I2.EQ.I3OLD .AND. I3.EQ.I1OLD
     2.OR.I1.EQ.I3OLD .AND. I2.EQ.I1OLD .AND. I3.EQ.I2OLD)
     3                                           GOTO 500
C
C     LOAD PARTIAL DERIVATIVES
C
      DO 100 J=1,5
        PDL(J,1)= PD(I1,J)
        PDL(J,2)= PD(I2,J)
        PDL(J,3)= PD(I3,J)
  100 CONTINUE
C
C     LOAD COORDINATES
C
      X(1)= XD(I1)
      Y(1)= YD(I1)
      X(2)= XD(I2)
      Y(2)= YD(I2)
      X(3)= XD(I3)
      Y(3)= YD(I3)
      Z(1)= ZD(I1)
      Z(2)= ZD(I2)
      Z(3)= ZD(I3)
C
C     COMPUTE COEFFICIENTS ALONG SIDES
      CALL C1SIDE (X,Y,Z,PDL,2)
C
C     COMPUTE COEFFICIENTS FOR INSIDE OF TRIANGLE
C     AND CONSTANTS FOR TRANSF. CARTESIAN TO TRIANGULAR
      CALL C1INSD (PDL)
C
C     EVALUATE POLYNOMIAL AT XP,YP
  500 CONTINUE
      XR= XP -X(3)
      YR= YP -Y(3)
      CALL C1HORN (XR,YR,ZP)
C
      IOUTS= 0
      GOTO 2000
C
C     POINT IS OUTSIDE TRIANGULATION
C
 1000 IOUTS= 1
      ZP= 0.
      IF (I3.NE.0) IST= I3
      IF (I2.NE.0) IST= I2
      IF (I1.NE.0) IST= I1
 3000 RETURN
C
C     SAVE VALUES OF I1,I2,I3
 2000 CONTINUE
      I1OLD= I1
      I2OLD= I2
      I3OLD= I3
C
      RETURN
      END
      SUBROUTINE C1SIDE (X,Y,Z,PDL,NSIDES)
C
C     C1-GRID
C
C     COMPUTATION OF COEFFICIENTS FOR POLYNOMIALS ALONG SIDES
C
C     AUTHOR      : A. PREUSSER, 1989
C                   FRITZ-HABER-INSTITUT
C                   DER MAX-PLANCK-GESELLSCHAFT
C                   FARADAYWEG 4-6
C                   D-1000 BERLIN 33
C
C
C     X        X-COORDINATES OF VERTICES
C     Y        Y-COORDINATES OF VERTICES
C     PDL(J,K) PARTIAL DERIVATIVES AT VERTEX K
C              (ZX,ZY,ZXX,ZXY,ZYY)
C     NSIDES   NUMBER OF POLYNOMIALS TO BE DETERMINED
C              2, FOR SIDE 1 AND 2
C              3, FOR SIDE 1, 2, AND 3
C
C                                  2
C                                  *
C                                /  .
C                              /     .
C                            /        .
C                          /           .
C                        /              .
C             SIDE(1)  /                 .
C                    /                    . SIDE(3)
C                  /                       .
C                /                          .
C              /                             .
C            /                                .
C          /                                   .
C        /                                      .
C    3 *-----------------------------------------* 1
C                      SIDE(2)
C
C
C
      COMMON /C1PCO/ P0(3),P1(3),P2(3),P3(3),P4(3),P5(3)
     1,              P11,P12,P13,P14,P21,P22,P23,P31,P32,P41
C     P0...P5     COEFFICIENTS OF POLYNOMIALS ALONG THE SIDES
C     P11...P41   COEFFICIENTS FOR BIVARIATE POLYNOMIAL INSIDE TRIANGLE
C
      COMMON /C1CCO/ DX(3),DY(3),DX2(3),DY2(3),SL2(3)
     1,              AP,BP,CP,DP
C
      COMMON /C1DCO/ ZS(3,3),ZSS(3,3)
C
C
C     OUTPUT OF THIS ROUTINE
C     ----------------------
C     DX,DY     COORDINATE DIFFERENCES OF ENDPOINTS OF SIDES
C     DX2,DY2   DX**2,DY**2
C     SL2       SIDE LENGTH **2
C     ZS        PARTIAL DERIVATIVES WITH RESPECT TO TRIANGULAR
C               COORDINATES L1,L2,L3 (CORRESPOND TO DIRECTIONS
C               OF SIDES 2,1,3 RESPECTIVELY)
C               FIRST INDEX: VARIABLE
C               SECOND INDEX: POINT
C     ZSS       SECOND PARTIAL DERIVATIVES
C               SIMILAR TO ZS
C     P0...P5   COEFFICIENTS FOR POLYNOMIALS ALONG THE THREE
C               SIDES
C
C
      DIMENSION PDL(5,3),X(3),Y(3),Z(3)
C
C     COMPUTE COORDINATE DIFFERENCES FOR SIDES
      DO 60 I=1,3
        NP1= 3
        NP2= 2
        IF (I.EQ.3) NP1= 1
        IF (I.EQ.2) NP2= 1
        DX(I)= X(NP2) - X(NP1)
        DY(I)= Y(NP2) - Y(NP1)
        DX2(I)= DX(I)*DX(I)
        DY2(I)= DY(I)*DY(I)
        SL2(I)= DX2(I) + DY2(I)
   60 CONTINUE
C
C     CONVERTS THE PARTIAL DERIVATIVES AT THE VERTICES OF THE
C     TRIANGLE FROM THE CARTESIAN X-Y-SYSTEM TO THE TRIANGULAR
C     L1-L2-SYSTEM USING THE CHAIN RULE
      DXDY1= 2.0*DX(1)*DY(1)
      DXDY2= 2.0*DX(2)*DY(2)
      DXDY3= 2.0*DX(3)*DY(3)
      DO 260 K=1,3
        ZS(1,K) = DX(2)*PDL(1,K)  + DY(2)*PDL(2,K)
        ZS(2,K) = DX(1)*PDL(1,K)  + DY(1)*PDL(2,K)
        ZSS(1,K)= DX2(2)*PDL(3,K) + DXDY2*PDL(4,K) + DY2(2)*PDL(5,K)
        ZSS(2,K)= DX2(1)*PDL(3,K) + DXDY1*PDL(4,K) + DY2(1)*PDL(5,K)
  260 CONTINUE
C
C     OPTIONALLY, COMPUTE PARTIAL DERIVATIVES FOR THIRD
C     VARIABLE L3 (AT POINTS 1 AND 2)
      IF (NSIDES.LT.3) GOTO 280
      DO 270 K=1,2
        ZS(3,K) = DX(3)*PDL(1,K)  + DY(3)*PDL(2,K)
        ZSS(3,K)= DX2(3)*PDL(3,K) + DXDY3*PDL(4,K) + DY2(3)*PDL(5,K)
  270 CONTINUE
  280 CONTINUE
C
C     CALCULATES THE COEFFICIENTS OF THE POLYNOMIALS ALONG
C     THE TWO (OR THREE) SIDES OF THE TRIANGLE
      DO 300 I=1,NSIDES
        NP1= 3
        NP2= 2
        IF (I.EQ.3) NP1= 1
        IF (I.EQ.2) NP2= 1
        J= NP2
        IF (I.EQ.3) J= 3
        P0(I)= Z(NP1)
        P1(I)= ZS(J,NP1)
        P2(I)= 0.5*ZSS(J,NP1)
        H1= Z(NP2) - P0(I) - P1(I) - P2(I)
        H2= ZS(J,NP2) - P1(I) - ZSS(J,NP1)
        H3= ZSS(J,NP2) - ZSS(J,NP1)
        P3(I)= 10.0*H1-4.0*H2+0.5*H3
        P4(I)=-15.0*H1+7.0*H2    -H3
        P5(I)=  6.0*H1-3.0*H2+0.5*H3
  300 CONTINUE
C
      RETURN
      END
      SUBROUTINE C1INSD (PDL)
C
C     C1-GRID
C
C     COMPUTATION OF POLYNOMIAL COEFFICIENTS P11...P41 FOR
C     BIVARIATE POLYNOMIAL INSIDE TRIANGLE
C     AND
C     CONSTANTS FOR TRANSFORMATION OF
C     CARTESIAN TO TRIANGULAR COORDINATES
C
C     AUTHOR      : A. PREUSSER, 1989
C                   FRITZ-HABER-INSTITUT
C                   DER MAX-PLANCK-GESELLSCHAFT
C                   FARADAYWEG 4-6
C                   D-1000 BERLIN 33
C
      DIMENSION PDL(5,3),Z12(3)
C
C     PDL     PARTIAL DERIVATIVES AT VERTICES 1,2,3
C             FIRST INDEX: 1, ZX
C                        : 2, ZY
C                        : 3, ZXX
C                        : 4, ZXY
C                        : 5, ZYY
C
C     Z12     MIXED PARTIAL DERIVATIVES WITH RESPECT TO L1,L2
C             AT THE THREE VERTICES 1,2,3
C
      COMMON /C1PCO/  P0(3),P1(3),P2(3),P3(3),P4(3),P5(3)
     3,     P11,P12,P13,P14,P21,P22,P23,P31,P32,P41
C     P0...P5     COEFFICIENTS OF POLYNOMIALS ALONG THE SIDES
C     P11...P41   COEFFICIENTS FOR BIVARIATE POLYNOMIAL INSIDE TRIANGLE
C
C
      COMMON /C1CCO/ DX(3),DY(3),DX2(3),DY2(3),SL2(3)
     1,              AP,BP,CP,DP
C
      COMMON /C1DCO/ ZS(3,3),ZSS(3,3)
C
C
C     INPUT TO THIS ROUTINE
C     ---------------------
C     VARIABLES FROM /C1CCO/ AND /C1DCO/ COMPUTED IN C1SIDE
C
C
C     OUTPUT OF THIS ROUTINE
C     ----------------------
C     P11...P41 POLYNOMIAL COEFFICIENTS FOR INSIDE REPRESENTATION
C     AP,BP,CP,DP    CONST. FOR TRANSF. CART. TO TRIANG. SYSTEM
C
C
C--------------------------------------------------------------------
C
C     COMPUTE CONST. FOR TRANSF. CART.-TRIANG. COORD. SYSTEM
      AD=  DX(2)*DY(1)
      BC=  DX(1)*DY(2)
      DLT= AD-BC
      AP=  DY(1)/DLT
      BP= -DX(1)/DLT
      CP= -DY(2)/DLT
      DP=  DX(2)/DLT
C
C     COMPUTE MIXED PARTIALS
      ADBC= AD+BC
      AB= DX(2)*DX(1)
      CD= DY(1)*DY(2)
      DO 100 K=1,3
  100   Z12(K)= AB*PDL(3,K) + ADBC*PDL(4,K) + CD*PDL(5,K)
C
C     COMPUTE COEFFICIENTS FOR REPRESENTATION INSIDE
C             TRIANGLE
      ABCD5= 5.*(AB+CD)
      P14= P5(1) * ABCD5/SL2(1)
      P41= P5(2) * ABCD5/SL2(2)
      P11= Z12(3)
      H1= ZS(2,1)-P1(1)-P11-P41
      H2= Z12(1)-P11-4.0*P41
      P21= 3.0*H1-H2
      P31=-2.0*H1+H2
      H1= ZS(1,2)-P1(2)-P11-P14
      H2= Z12(2)-P11-4.0*P14
      P12= 3.0*H1-H2
      P13=-2.0*H1+H2
      H1= 0.5*ZSS(2,1)-P2(1)-P12
      H2= 0.5*ZSS(1,2)-P2(2)-P21
      E1= 2.5*(SL2(2)-SL2(1))/SL2(3) + 2.5
C     EQUIVALENT TO
C     E1= -5.*(DX(2)*DX(3)+DY(2)*DY(3))/SL2(3)
      G1= 3. - E1
      G2=-2. + E1
      G3= E1*(P5(1) - P5(2) + P41 -P14)
     1    + P14 - 4.*P41 + 5.*P5(2)
      P22= G1*H1 + G2*H2 + G3
      P32= H1-P22
      P23= H2-P22
      RETURN
      END
      SUBROUTINE C1HORN (X,Y,Z)
C
C     C1-GRID
C
C     EVALUATION OF BIVARIATE QUINTIC POLYNOMIAL
C
C     AUTHOR      : A. PREUSSER, 1989
C                   FRITZ-HABER-INSTITUT
C                   DER MAX-PLANCK-GESELLSCHAFT
C                   FARADAYWEG 4-6
C                   D-1000 BERLIN 33
C
C
C    INPUT PARAMETERS
C    ----------------
C    X,Y      CARTESIAN COORDINATES RELATIVE TO VERTEX (3)
C             OF TRIANGLE
C
C    OUTPUT PARAMETER
C    ----------------
C    Z        VALUE OF POLYNOMIAL AT X,Y
C
C
      COMMON /C1PCO/  P0(3),P1(3),P2(3),P3(3),P4(3),P5(3)
     3,           P11,P12,P13,P14,P21,P22,P23,P31,P32,P41
C
      COMMON /C1CCO/ DX(3),DY(3),DX2(3),DY2(3),SL2(3)
     1,              AP,BP,CP,DP
C
C     INPUT TO THIS ROUTINE
C     ---------------------
C     ALL VARIABLES FROM /C1PCO/ (POLYNOMIAL COEFFICIENTS)
C     AP,BP,CP,DP    CONSTANTS DERIVED FROM COORDINATE DIFFERENCES
C                    ON SIDES
C
C ------------------------------------------------------------------
C
C     TRANSFORMATION FROM CARTESIAN TO TRIANGULAR COORDINATES
C                                                (U=L1, V=L2)
      U= AP*X + BP*Y
      V= CP*X + DP*Y
C
C     EVALUATION BY HORNER'S SCHEME (BIVARIATE)
      H0=P0(1)+V*(P1(1)+V*(P2(1)+V*(P3(1)+V*(P4(1)+V*P5(1)))))
      H1=P1(2)+V*(P11+V*(P12+V*(P13+V*P14)))
      H2=P2(2)+V*(P21+V*(P22+V*P23))
      H3=P3(2)+V*(P31+V*P32)
      H4=P4(2)+V*P41
      Z= H0+U*(H1+U*(H2+U*(H3+U*(H4+U*P5(2)))))
C
      RETURN
      END
      SUBROUTINE GRANKA(N,X,Y,Z,IADJ,IEND,EPS, NIT,
     .                  ZXZY, IER)
      INTEGER N, IADJ(*), IEND(N), NIT, IER
      REAL    X(N), Y(N), Z(N), EPS, ZXZY(N,2)
C
C***********************************************************
C
C                                               ROBERT RENKA
C                                       OAK RIDGE NATL. LAB.
C                                             (615) 576-5139
C
C *  THIS PROGRAM WAS MODIFIED BY A.PREUSSER, 1989.
C *  1)
C *  FIRST AND SECOND INDEX OF ZXZY ARE INTERCHANGED.
C *  NOTE THE REMARK BELOW.
C *  2)
C *  EPS IS CHANGED FROM AN ABSOLUTE TO A RELATIVE ERROR.
C
C   GIVEN A TRIANGULATION OF N NODES IN THE PLANE WITH
C ASSOCIATED DATA VALUES, THIS ROUTINE USES A GLOBAL METHOD
C TO COMPUTE ESTIMATED GRADIENTS AT THE NODES.  THE METHOD
C CONSISTS OF MINIMIZING A QUADRATIC FUNCTIONAL Q(G) OVER
C THE N-VECTOR G OF GRADIENTS WHERE Q APPROXIMATES THE LIN-
C EARIZED CURVATURE OF AN INTERPOLANT F OVER THE TRIANGULA-
C TION.  THE RESTRICTION OF F TO AN ARC OF THE TRIANGULATION
C IS TAKEN TO BE THE HERMITE CUBIC INTERPOLANT OF THE DATA
C VALUES AND TANGENTIAL GRADIENT COMPONENTS AT THE END-
C POINTS OF THE ARC, AND Q IS THE SUM OF THE LINEARIZED
C CURVATURES OF F ALONG THE ARCS -- THE INTEGRALS OVER THE
C ARCS OF D2F(T)**2 WHERE D2F(T) IS THE SECOND DERIVATIVE
C OF F WITH RESPECT TO DISTANCE T ALONG THE ARC.  THIS MIN-
C IMIZATION PROBLEM CORRESPONDS TO AN ORDER 2N SYMMETRIC
C POSITIVE-DEFINITE SPARSE LINEAR SYSTEM WHICH IS SOLVED FOR
C THE X AND Y PARTIAL DERIVATIVES BY THE BLOCK GAUSS-SEIDEL
C METHOD WITH 2 BY 2 BLOCKS.
C   AN ALTERNATIVE METHOD, SUBROUTINE GRADL, COMPUTES A
C LOCAL APPROXIMATION TO THE PARTIALS AT A SINGLE NODE AND
C MAY BE MORE ACCURATE, DEPENDING ON THE DATA VALUES AND
C DISTRIBUTION OF NODES (NEITHER METHOD EMERGED AS SUPERIOR
C IN TESTS FOR ACCURACY).  HOWEVER, IN TESTS RUN ON AN IBM
C 370, GRANKA WAS FOUND TO BE ABOUT 3.6 TIMES AS FAST FOR
C NIT = 4.
C
C INPUT PARAMETERS - N - NUMBER OF NODES.  N .GE. 3.
C
C                  X,Y - CARTESIAN COORDINATES OF THE NODES.
C
C                    Z - DATA VALUES AT THE NODES.  Z(I) IS
C                        ASSOCIATED WITH (X(I),Y(I)).
C
C            IADJ,IEND - DATA STRUCTURE DEFINING THE TRIAN-
C                        GULATION.  SEE SUBROUTINE TRMESH.
C
C                  EPS - NONNEGATIVE CONVERGENCE CRITERION.
C                        THE METHOD IS TERMINATED WHEN THE
C                        MAXIMUM CHANGE IN A GRADIENT COMPO-
C                        NENT BETWEEN ITERATIONS IS AT MOST
C                        EPS.  EPS = 1.E-2 IS SUFFICIENT FOR
C                        EFFECTIVE CONVERGENCE.
C *                      EPS IS RELATIVE TO THE UPDATED COMPONENT.
C
C THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
C
C                  NIT - MAXIMUM NUMBER OF GAUSS-SEIDEL
C                        ITERATIONS TO BE APPLIED.  THIS
C                        MAXIMUM WILL LIKELY BE ACHIEVED IF
C                        EPS IS SMALLER THAN THE MACHINE
C                        PRECISION.  OPTIMAL EFFICIENCY WAS
C                        ACHIEVED IN TESTING WITH EPS = 0
C                        AND NIT = 3 OR 4.
C
C *               ZXZY - N BY 2 ARRAY WHOSE COLUMNS CONTAIN
C                        INITIAL ESTIMATES OF THE PARTIAL
C                        DERIVATIVES (ZERO VECTORS ARE
C                        SUFFICIENT).
C
C *                      NOTE THAT ZXZY MUST BE DECLARED
C *                      EXACTLY (N,2) IN THE CALLING
C *                      PROGRAM IN ORDER TO ENSHURE
C *                      CORRECT ADDRESSING.
C
C OUTPUT PARAMETERS - NIT - NUMBER OF GAUSS-SEIDEL ITERA-
C                           TIONS EMPLOYED.
C
C                    ZXZY - ESTIMATED X AND Y PARTIAL DERIV-
C                           ATIVES AT THE NODES WITH X PAR-
C *                         TIALS IN THE FIRST COL.  ZXZY IS
C                           NOT CHANGED IF IER = 2.
C
C                     IER - ERROR INDICATOR
C                           IER = 0 IF THE CONVERGENCE CRI-
C                                   TERION WAS ACHIEVED.
C                           IER = 1 IF CONVERGENCE WAS NOT
C                                   ACHIEVED WITHIN NIT
C                                   ITERATIONS.
C                           IER = 2 IF N OR EPS IS OUT OF
C                                   RANGE OR NIT .LT. 0 ON
C                                   INPUT.
C
C MODULES REFERENCED BY GRANKA - NONE
C
C INTRINSIC FUNCTIONS CALLED BY GRANKA - SQRT, AMAX1, ABS
C
C***********************************************************
C
      INTEGER NN, MAXIT, ITER, K, INDF, INDL, INDX, NB
      REAL    TOL, DMAREL, XK, YK, ZK, ZXK, ZYK, A11, A12,
     .        A22, R1, R2, DELX, DELY, DELXS, DELYS, DSQ,
     .        DCUB, T, DZX, DZY
C
C LOCAL PARAMETERS -
C
C NN =          LOCAL COPY OF N
C MAXIT =       INPUT VALUE OF NIT
C ITER =        NUMBER OF ITERATIONS USED
C K =           DO-LOOP AND NODE INDEX
C INDF,INDL =   IADJ INDICES OF THE FIRST AND LAST NEIGHBORS
C                 OF K
C INDX =        IADJ INDEX IN THE RANGE INDF,...,INDL
C NB =          NEIGHBOR OF K
C TOL =         LOCAL COPY OF EPS
C DMAREL =      MAXIMUM CHANGE IN A GRADIENT COMPONENT BE-
C *               TWEEN ITERATIONS RELATIVE TO THE UPDATED
C *               COMPONENT
C XK,YK,ZK =    X(K), Y(K), Z(K)
C ZXK,ZYK =  *  INITIAL VALUES OF ZXZY(K,1) AND ZXZY(K,2)
C A11,A12,A22 = MATRIX COMPONENTS OF THE 2 BY 2 BLOCK A*DG
C                 = R WHERE A IS SYMMETRIC, DG = (DZX,DZY)
C                 IS THE CHANGE IN THE GRADIENT AT K, AND R
C                 IS THE RESIDUAL
C R1,R2 =       COMPONENTS OF THE RESIDUAL -- DERIVATIVES OF
C                 Q WITH RESPECT TO THE COMPONENTS OF THE
C                 GRADIENT AT NODE K
C DELX,DELY =   COMPONENTS OF THE ARC NB-K
C DELXS,DELYS = DELX**2, DELY**2
C DSQ =         SQUARE OF THE DISTANCE D BETWEEN K AND NB
C DCUB =        D**3
C T =           FACTOR OF R1 AND R2
C DZX,DZY =     SOLUTION OF THE 2 BY 2 SYSTEM -- CHANGE IN
C                 DERIVATIVES AT K FROM THE PREVIOUS ITERATE
C
      NN = N
      TOL = EPS
      MAXIT = NIT
C
C ERROR CHECKS AND INITIALIZATION
C
      IF (NN .LT. 3  .OR.  TOL .LT. 0.  .OR.  MAXIT .LT. 0)
     .   GO TO 5
      ITER = 0
C
C TOP OF ITERATION LOOP
C
    1 IF (ITER .EQ. MAXIT) GO TO 4
      DMAREL = 0.
      INDL = 0
      DO 3 K = 1,NN
        XK = X(K)
        YK = Y(K)
        ZK = Z(K)
        ZXK = ZXZY(K,1)
        ZYK = ZXZY(K,2)
C
C   INITIALIZE COMPONENTS OF THE 2 BY 2 SYSTEM
C
        A11 = 0.
        A12 = 0.
        A22 = 0.
        R1 = 0.
        R2 = 0.
C
C   LOOP ON NEIGHBORS NB OF K
C
        INDF = INDL + 1
        INDL = IEND(K)
        DO 2 INDX = INDF,INDL
          NB = IADJ(INDX)
          IF (NB .EQ. 0) GO TO 2
C
C   COMPUTE THE COMPONENTS OF ARC NB-K
C
          DELX = X(NB) - XK
          DELY = Y(NB) - YK
          DELXS = DELX*DELX
          DELYS = DELY*DELY
          DSQ = DELXS + DELYS
          DCUB = DSQ*SQRT(DSQ)
C
C   UPDATE THE SYSTEM COMPONENTS FOR NODE NB
C
          A11 = A11 + DELXS/DCUB
          A12 = A12 + DELX*DELY/DCUB
          A22 = A22 + DELYS/DCUB
          T = ( 1.5*(Z(NB)-ZK) - ((ZXZY(NB,1)/2.+ZXK)*DELX +
     .          (ZXZY(NB,2)/2.+ZYK)*DELY) )/DCUB
          R1 = R1 + T*DELX
          R2 = R2 + T*DELY
    2     CONTINUE
C
C * SOLVE THE 2 BY 2 SYSTEM
C
        DZY = (A11*R2 - A12*R1)/(A11*A22 - A12*A12)
        DZX = (R1 - A12*DZY)/A11
C
C   UPDATE THE PARTIALS AT NODE K
C
        ZXZY(K,1) = ZXK + DZX
        ZXZY(K,2) = ZYK + DZY
C
C * UPDATE DMAREL
C
        ZXREL= 1.
        ZYREL= 1.
        IF (ZXK.NE.0.)  ZXREL= ABS(DZX/ZXK)
        IF (ZYK.NE.0.)  ZYREL= ABS(DZY/ZYK)
        DMAREL= AMAX1(DMAREL,ZXREL,ZYREL)
C
    3 CONTINUE
*
C
C   INCREMENT ITER AND TEST FOR CONVERGENCE
C
      ITER = ITER + 1
      IF (DMAREL .GT. TOL) GO TO 1
C
C METHOD CONVERGED
C
      NIT = ITER
      IER = 0
      RETURN
C
C METHOD FAILED TO CONVERGE WITHIN NIT ITERATIONS
C
    4 IER = 1
      RETURN
C
C PARAMETER OUT OF RANGE
C
    5 NIT = 0
      IER = 2
      RETURN
      END
      SUBROUTINE TRMESH (N,X,Y, IADJ,IEND,IER)
      INTEGER N, IADJ(*), IEND(N), IER
      REAL    X(N), Y(N)
C
C***********************************************************
C
C                                               ROBERT RENKA
C                                       OAK RIDGE NATL. LAB.
C                                             (615) 576-5139
C
C   THIS ROUTINE CREATES A THIESSEN TRIANGULATION OF N
C ARBITRARILY SPACED POINTS IN THE PLANE REFERRED TO AS
C NODES.  THE TRIANGULATION IS OPTIMAL IN THE SENSE THAT IT
C IS AS NEARLY EQUIANGULAR AS POSSIBLE.  TRMESH IS PART OF
C AN INTERPOLATION PACKAGE WHICH ALSO PROVIDES SUBROUTINES
C TO REORDER THE NODES, ADD A NEW NODE, DELETE AN ARC, PLOT
C THE MESH, AND PRINT THE DATA STRUCTURE.
C   UNLESS THE NODES ARE ALREADY ORDERED IN SOME REASONABLE
C FASHION, THEY SHOULD BE REORDERED BY SUBROUTINE REORDR FOR
C INCREASED EFFICIENCY BEFORE CALLING TRMESH.
C
C INPUT PARAMETERS -     N - NUMBER OF NODES IN THE MESH.
C                            N .GE. 3.
C
C                      X,Y - N-VECTORS OF COORDINATES.
C                            (X(I),Y(I)) DEFINES NODE I.
C
C                     IADJ - VECTOR OF LENGTH .GE. 6*N-9.
C
C                     IEND - VECTOR OF LENGTH .GE. N.
C
C N, X, AND Y ARE NOT ALTERED BY THIS ROUTINE.
C
C OUTPUT PARAMETERS - IADJ - ADJACENCY LISTS OF NEIGHBORS IN
C                            COUNTERCLOCKWISE ORDER.  THE
C                            LIST FOR NODE I+1 FOLLOWS THAT
C                            FOR NODE I WHERE X AND Y DEFINE
C                            THE ORDER.  THE VALUE 0 DENOTES
C                            THE BOUNDARY (OR A PSEUDO-NODE
C                            AT INFINITY) AND IS ALWAYS THE
C                            LAST NEIGHBOR OF A BOUNDARY
C                            NODE.  IADJ IS UNCHANGED IF IER
C                            .NE. 0.
C
C                     IEND - POINTERS TO THE ENDS OF
C                            ADJACENCY LISTS (SETS OF
C                            NEIGHBORS) IN IADJ.  THE
C                            NEIGHBORS OF NODE 1 BEGIN IN
C                            IADJ(1).  FOR K .GT. 1, THE
C                            NEIGHBORS OF NODE K BEGIN IN
C                            IADJ(IEND(K-1)+1) AND K HAS
C                            IEND(K) - IEND(K-1) NEIGHBORS
C                            INCLUDING (POSSIBLY) THE
C                            BOUNDARY.  IADJ(IEND(K)) .EQ. 0
C                            IFF NODE K IS ON THE BOUNDARY.
C                            IEND IS UNCHANGED IF IER = 1.
C                            IF IER = 2 IEND CONTAINS THE
C                            INDICES OF A SEQUENCE OF N
C                            NODES ORDERED FROM LEFT TO
C                            RIGHT WHERE LEFT AND RIGHT ARE
C                            DEFINED BY ASSUMING NODE 1 IS
C                            TO THE LEFT OF NODE 2.
C
C                      IER - ERROR INDICATOR
C                            IER = 0 IF NO ERRORS WERE
C                                    ENCOUNTERED.
C                            IER = 1 IF N .LT. 3.
C                            IER = 2 IF N .GE. 3 AND ALL
C                                    NODES ARE COLLINEAR.
C
C MODULES REFERENCED BY TRMESH - SHIFTD, ADNODE, TRFIND,
C                                INTADD, BDYADD, SWPTST,
C                                SWAP, INDEX
C
C***********************************************************
C
      INTEGER NN, K, KM1, NL, NR, IND, INDX, N0, ITEMP,
     .        IERR, KM1D2, KMI, I, KMIN
      REAL    XL, YL, XR, YR, DXR, DYR, XK, YK, DXK, DYK,
     .        CPROD, SPROD
C
C LOCAL PARAMETERS -
C
C NN =          LOCAL COPY OF N
C K =           NODE (INDEX) TO BE INSERTED INTO IEND
C KM1 =         K-1 - (VARIABLE) LENGTH OF IEND
C NL,NR =       IEND(1), IEND(KM1) -- LEFTMOST AND RIGHTMOST
C                 NODES IN IEND AS VIEWED FROM THE RIGHT OF
C                 1-2 WHEN IEND CONTAINS THE INITIAL ORDERED
C                 SET OF NODAL INDICES
C XL,YL,XR,YR = X AND Y COORDINATES OF NL AND NR
C DXR,DYR =     XR-XL, YR-YL
C XK,YK =       X AND Y COORDINATES OF NODE K
C DXK,DYK =     XK-XL, YK-YL
C CPROD =       VECTOR CROSS PRODUCT OF NL-NR AND NL-K --
C                 USED TO DETERMINE THE POSITION OF NODE K
C                 WITH RESPECT TO THE LINE DEFINED BY THE
C                 NODES IN IEND
C SPROD =       SCALAR PRODUCT USED TO DETERMINE THE
C                 INTERVAL CONTAINING NODE K WHEN K IS ON
C                 THE LINE DEFINED BY THE NODES IN IEND
C IND,INDX =    INDICES FOR IEND AND IADJ, RESPECTIVELY
C N0,ITEMP =    TEMPORARY NODES (INDICES)
C IERR =        DUMMY PARAMETER FOR CALL TO ADNODE
C KM1D2,KMI,I = KM1/2, K-I, DO-LOOP INDEX -- USED IN IEND
C                 REORDERING LOOP
C KMIN =        FIRST NODE INDEX SENT TO ADNODE
C
      NN = N
      IER = 1
      IF (NN .LT. 3) RETURN
      IER = 0
C
C INITIALIZE IEND, NL, NR, AND K
C
      IEND(1) = 1
      IEND(2) = 2
      XL = X(1)
      YL = Y(1)
      XR = X(2)
      YR = Y(2)
      K = 2
C
C BEGIN LOOP ON NODES 3,4,...
C
    1 DXR = XR-XL
      DYR = YR-YL
C
C NEXT LOOP BEGINS HERE IF NL AND NR ARE UNCHANGED
C
    2 IF (K .EQ. NN) GO TO 13
      KM1 = K
      K = KM1 + 1
      XK = X(K)
      YK = Y(K)
      DXK = XK-XL
      DYK = YK-YL
      CPROD = DXR*DYK - DXK*DYR
      IF (CPROD .GT. 0.) GO TO 6
      IF (CPROD .LT. 0.) GO TO 8
C
C NODE K LIES ON THE LINE CONTAINING NODES 1,2,...,K-1.
C   SET SPROD TO (NL-NR,NL-K).
C
      SPROD = DXR*DXK + DYR*DYK
      IF (SPROD .GT. 0.) GO TO 3
C
C NODE K IS TO THE LEFT OF NL.  INSERT K AS THE FIRST
C   (LEFTMOST) NODE IN IEND AND SET NL TO K.
C
      CALL SHIFTD(1,KM1,1, IEND )
      IEND(1) = K
      XL = XK
      YL = YK
      GO TO 1
C
C NODE K IS TO THE RIGHT OF NL.  FIND THE LEFTMOST NODE
C   N0 WHICH LIES TO THE RIGHT OF K.
C   SET SPROD TO (N0-NL,N0-K).
C
    3 DO 4 IND = 2,KM1
        N0 = IEND(IND)
        SPROD = (XL-X(N0))*(XK-X(N0)) +
     .          (YL-Y(N0))*(YK-Y(N0))
        IF (SPROD .GE. 0.) GO TO 5
    4   CONTINUE
C
C NODE K IS TO THE RIGHT OF NR.  INSERT K AS THE LAST
C   (RIGHTMOST) NODE IN IEND AND SET NR TO K.
C
      IEND(K) = K
      XR = XK
      YR = YK
      GO TO 1
C
C NODE K LIES BETWEEN IEND(IND-1) AND IEND(IND).  INSERT K
C   IN IEND.
C
    5 CALL SHIFTD(IND,KM1,1, IEND )
      IEND(IND) = K
      GO TO 2
C
C NODE K IS TO THE LEFT OF NL-NR.  REORDER IEND SO THAT NL
C   IS THE LEFTMOST NODE AS VIEWED FROM K.
C
    6 KM1D2 = KM1/2
      DO 7 I = 1,KM1D2
        KMI = K-I
        ITEMP = IEND(I)
        IEND(I) = IEND(KMI)
        IEND(KMI) = ITEMP
    7   CONTINUE
C
C NODE K IS TO THE RIGHT OF NL-NR.  CREATE A TRIANGULATION
C   CONSISTING OF NODES 1,2,...,K.
C
    8 NL = IEND(1)
      NR = IEND(KM1)
C
C CREATE THE ADJACENCY LISTS FOR THE FIRST K-1 NODES.
C   INSERT NEIGHBORS IN REVERSE ORDER.  EACH NODE HAS FOUR
C   NEIGHBORS EXCEPT NL AND NR WHICH HAVE THREE.
C
      DO 9 IND = 1,KM1
        N0 = IEND(IND)
        INDX = 4*N0
        IF (N0 .GE. NL) INDX = INDX-1
        IF (N0 .GE. NR) INDX = INDX-1
        IADJ(INDX) = 0
        INDX = INDX-1
        IF (IND .LT. KM1) IADJ(INDX) = IEND(IND+1)
        IF (IND .LT. KM1) INDX = INDX-1
        IADJ(INDX) = K
        IF (IND .EQ. 1) GO TO 9
        IADJ(INDX-1) = IEND(IND-1)
    9   CONTINUE
C
C CREATE THE ADJACENCY LIST FOR NODE K
C
      INDX = 5*KM1 - 1
      IADJ(INDX) = 0
      DO 10 IND = 1,KM1
        INDX = INDX-1
        IADJ(INDX) = IEND(IND)
   10   CONTINUE
C
C REPLACE IEND ELEMENTS WITH POINTERS TO IADJ
C
      INDX = 0
      DO 11 IND = 1,KM1
        INDX = INDX + 4
        IF (IND .EQ. NL  .OR.  IND .EQ. NR) INDX = INDX-1
        IEND(IND) = INDX
   11   CONTINUE
      INDX = INDX + K
      IEND(K) = INDX
C
C ADD THE REMAINING NODES TO THE TRIANGULATION
C
      IF (K .EQ. NN) RETURN
      KMIN = K+1
      DO 12 K = KMIN,NN
        CALL ADNODE(K,X,Y, IADJ,IEND, IERR)
   12   CONTINUE
      RETURN
C
C ALL NODES ARE COLLINEAR
C
   13 IER = 2
      RETURN
      END
      SUBROUTINE ADNODE (KK,X,Y, IADJ,IEND, IER)
      INTEGER KK, IADJ(*), IEND(KK), IER
      REAL    X(KK), Y(KK)
      LOGICAL SWPTST
      EXTERNAL INDEX
C
C***********************************************************
C
C                                               ROBERT RENKA
C                                       OAK RIDGE NATL. LAB.
C                                             (615) 576-5139
C
C   THIS ROUTINE ADDS NODE KK TO A TRIANGULATION OF A SET
C OF POINTS IN THE PLANE PRODUCING A NEW TRIANGULATION.  A
C SEQUENCE OF EDGE SWAPS IS THEN APPLIED TO THE MESH,
C RESULTING IN AN OPTIMAL TRIANGULATION.  ADNODE IS PART
C OF AN INTERPOLATION PACKAGE WHICH ALSO PROVIDES ROUTINES
C TO INITIALIZE THE DATA STRUCTURE, PLOT THE MESH, AND
C DELETE ARCS.
C
C INPUT PARAMETERS -   KK - INDEX OF THE NODE TO BE ADDED
C                           TO THE MESH.  KK .GE. 4.
C
C                     X,Y - VECTORS OF COORDINATES OF THE
C                           NODES IN THE MESH.  (X(I),Y(I))
C                           DEFINES NODE I FOR I = 1,..,KK.
C
C                    IADJ - SET OF ADJACENCY LISTS OF NODES
C                           1,..,KK-1.
C
C                    IEND - POINTERS TO THE ENDS OF
C                           ADJACENCY LISTS IN IADJ FOR
C                           EACH NODE IN THE MESH.
C
C IADJ AND IEND MAY BE CREATED BY TRMESH.
C
C KK, X, AND Y ARE NOT ALTERED BY THIS ROUTINE.
C
C OUTPUT PARAMETERS - IADJ,IEND - UPDATED WITH THE ADDITION
C                                 OF NODE KK AS THE LAST
C                                 ENTRY.
C
C                           IER - ERROR INDICATOR
C                                 IER = 0 IF NO ERRORS
C                                         WERE ENCOUNTERED.
C                                 IER = 1 IF ALL NODES
C                                         (INCLUDING KK) ARE
C                                         COLLINEAR.
C
C MODULES REFERENCED BY ADNODE - TRFIND, INTADD, BDYADD,
C                                SHIFTD, INDEX, SWPTST,
C                                SWAP
C
C***********************************************************
C
      INTEGER K, KM1, I1, I2, I3, INDKF, INDKL, NABOR1,
     .        IO1, IO2, IN1, INDK1, IND2F, IND21
      REAL    XK, YK
C
C LOCAL PARAMETERS -
C
C K =        LOCAL COPY OF KK
C KM1 =      K - 1
C I1,I2,I3 = VERTICES OF A TRIANGLE CONTAINING K
C INDKF =    IADJ INDEX OF THE FIRST NEIGHBOR OF K
C INDKL =    IADJ INDEX OF THE LAST NEIGHBOR OF K
C NABOR1 =   FIRST NEIGHBOR OF K BEFORE ANY SWAPS OCCUR
C IO1,IO2 =  ADJACENT NEIGHBORS OF K DEFINING AN ARC TO
C              BE TESTED FOR A SWAP
C IN1 =      VERTEX OPPOSITE K -- FIRST NEIGHBOR OF IO2
C              WHICH PRECEDES IO1.  IN1,IO1,IO2 ARE IN
C              COUNTERCLOCKWISE ORDER.
C INDK1 =    INDEX OF IO1 IN THE ADJACENCY LIST FOR K
C IND2F =    INDEX OF THE FIRST NEIGHBOR OF IO2
C IND21 =    INDEX OF IO1 IN THE ADJACENCY LIST FOR IO2
C XK,YK =    X(K), Y(K)
C
      IER = 0
      K = KK
C
C INITIALIZATION
C
      KM1 = K - 1
      XK = X(K)
      YK = Y(K)
C
C ADD NODE K TO THE MESH
C
      CALL TRFIND(KM1,XK,YK,X,Y,IADJ,IEND, I1,I2,I3)
      IF (I1 .EQ. 0) GO TO 5
      IF (I3 .EQ. 0) CALL BDYADD(K,I1,I2, IADJ,IEND )
      IF (I3 .NE. 0) CALL INTADD(K,I1,I2,I3, IADJ,IEND )
C
C INITIALIZE VARIABLES FOR OPTIMIZATION OF THE MESH
C
      INDKF = IEND(KM1) + 1
      INDKL = IEND(K)
      NABOR1 = IADJ(INDKF)
      IO2 = NABOR1
      INDK1 = INDKF + 1
      IO1 = IADJ(INDK1)
C
C BEGIN LOOP -- FIND THE VERTEX OPPOSITE K
C
    1 IND2F = 1
      IF (IO2 .NE. 1) IND2F = IEND(IO2-1) + 1
      IND21 = INDEX(IO2,IO1,IADJ,IEND)
      IF (IND2F .EQ. IND21) GO TO 2
      IN1 = IADJ(IND21-1)
      GO TO 3
C
C IN1 IS THE LAST NEIGHBOR OF IO2
C
    2 IND21 = IEND(IO2)
      IN1 = IADJ(IND21)
      IF (IN1 .EQ. 0) GO TO 4
C
C SWAP TEST -- IF A SWAP OCCURS, TWO NEW ARCS ARE OPPOSITE K
C              AND MUST BE TESTED.  INDK1 AND INDKF MUST BE
C              DECREMENTED.
C
    3 IF ( .NOT. SWPTST(IN1,K,IO1,IO2,X,Y) ) GO TO 4
      CALL SWAP(IN1,K,IO1,IO2, IADJ,IEND )
      IO1 = IN1
      INDK1 = INDK1 - 1
      INDKF = INDKF - 1
      GO TO 1
C
C NO SWAP OCCURRED.  RESET IO2 AND IO1, AND TEST FOR
C   TERMINATION.
C
    4 IF (IO1 .EQ. NABOR1) RETURN
      IO2 = IO1
      INDK1 = INDK1 + 1
      IF (INDK1 .GT. INDKL) INDK1 = INDKF
      IO1 = IADJ(INDK1)
      IF (IO1 .NE. 0) GO TO 1
      RETURN
C
C ALL NODES ARE COLLINEAR
C
    5 IER = 1
      RETURN
      END
      SUBROUTINE SHIFTD (NFRST,NLAST,KK, IARR )
      INTEGER NFRST, NLAST, KK, IARR(*)
C
C***********************************************************
C
C                                               ROBERT RENKA
C                                       OAK RIDGE NATL. LAB.
C                                             (615) 576-5139
C
C   THIS ROUTINE SHIFTS A SET OF CONTIGUOUS ELEMENTS OF AN
C INTEGER ARRAY KK POSITIONS DOWNWARD (UPWARD IF KK .LT. 0).
C THE LOOPS ARE UNROLLED IN ORDER TO INCREASE EFFICIENCY.
C
C INPUT PARAMETERS - NFRST,NLAST - BOUNDS ON THE PORTION OF
C                                  IARR TO BE SHIFTED.  ALL
C                                  ELEMENTS BETWEEN AND
C                                  INCLUDING THE BOUNDS ARE
C                                  SHIFTED UNLESS NFRST .GT.
C                                  NLAST, IN WHICH CASE NO
C                                  SHIFT OCCURS.
C
C                             KK - NUMBER OF POSITIONS EACH
C                                  ELEMENT IS TO BE SHIFTED.
C                                  IF KK .LT. 0 SHIFT UP.
C                                  IF KK .GT. 0 SHIFT DOWN.
C
C                           IARR - INTEGER ARRAY OF LENGTH
C                                  .GE. NLAST + MAX(KK,0).
C
C NFRST, NLAST, AND KK ARE NOT ALTERED BY THIS ROUTINE.
C
C OUTPUT PARAMETER -        IARR - SHIFTED ARRAY.
C
C MODULES REFERENCED BY SHIFTD - NONE
C
C***********************************************************
C
      INTEGER INC, K, NF, NL, NLP1, NS, NSL, I, IBAK, INDX,
     .        IMAX
      DATA    INC/5/
C
C LOCAL PARAMETERS -
C
C INC =  DO-LOOP INCREMENT (UNROLLING FACTOR) -- IF INC IS
C          CHANGED, STATEMENTS MUST BE ADDED TO OR DELETED
C          FROM THE DO-LOOPS
C K =    LOCAL COPY OF KK
C NF =   LOCAL COPY OF NFRST
C NL =   LOCAL COPY OF NLAST
C NLP1 = NL + 1
C NS =   NUMBER OF SHIFTS
C NSL =  NUMBER OF SHIFTS DONE IN UNROLLED DO-LOOP (MULTIPLE
C          OF INC)
C I =    DO-LOOP INDEX AND INDEX FOR IARR
C IBAK = INDEX FOR DOWNWARD SHIFT OF IARR
C INDX = INDEX FOR IARR
C IMAX = BOUND ON DO-LOOP INDEX
C
      K = KK
      NF = NFRST
      NL = NLAST
      IF (NF .GT. NL  .OR.  K .EQ. 0) RETURN
      NLP1 = NL + 1
      NS = NLP1 - NF
      NSL = INC*(NS/INC)
      IF ( K .LT. 0) GO TO 4
C
C SHIFT DOWNWARD STARTING FROM THE BOTTOM
C
      IF (NSL .LE. 0) GO TO 2
      DO 1 I = 1,NSL,INC
        IBAK = NLP1 - I
        INDX = IBAK + K
        IARR(INDX) = IARR(IBAK)
        IARR(INDX-1) = IARR(IBAK-1)
        IARR(INDX-2) = IARR(IBAK-2)
        IARR(INDX-3) = IARR(IBAK-3)
        IARR(INDX-4) = IARR(IBAK-4)
    1   CONTINUE
C
C PERFORM THE REMAINING NS-NSL SHIFTS ONE AT A TIME
C
    2 IBAK = NLP1 - NSL
    3 IF (IBAK .LE. NF) RETURN
      IBAK = IBAK - 1
      INDX = IBAK + K
      IARR(INDX) = IARR(IBAK)
      GO TO 3
C
C SHIFT UPWARD STARTING FROM THE TOP
C
    4 IF (NSL .LE. 0) GO TO 6
      IMAX = NLP1 - INC
      DO 5 I = NF,IMAX,INC
        INDX = I + K
        IARR(INDX) = IARR(I)
        IARR(INDX+1) = IARR(I+1)
        IARR(INDX+2) = IARR(I+2)
        IARR(INDX+3) = IARR(I+3)
        IARR(INDX+4) = IARR(I+4)
    5   CONTINUE
C
C PERFORM THE REMAINING NS-NSL SHIFTS ONE AT A TIME
C
    6 I = NSL + NF
    7 IF (I .GT. NL) RETURN
      INDX = I + K
      IARR(INDX) = IARR(I)
      I = I + 1
      GO TO 7
      END
      SUBROUTINE BDYADD (KK,I1,I2, IADJ,IEND )
      INTEGER KK, I1, I2, IADJ(*), IEND(KK)
C
C***********************************************************
C
C                                               ROBERT RENKA
C                                       OAK RIDGE NATL. LAB.
C                                             (615) 576-5139
C
C   THIS ROUTINE ADDS A BOUNDARY NODE TO A TRIANGULATION OF
C A SET OF KK-1 POINTS ON THE UNIT SPHERE.  IADJ AND IEND
C ARE UPDATED WITH THE INSERTION OF NODE KK.
C
C INPUT PARAMETERS -   KK - INDEX OF AN EXTERIOR NODE TO BE
C                           ADDED.  KK .GE. 4.
C
C                      I1 - FIRST (RIGHTMOST AS VIEWED FROM
C                           KK) BOUNDARY NODE IN THE MESH
C                           WHICH IS VISIBLE FROM KK - THE
C                           LINE SEGMENT KK-I1 INTERSECTS
C                           NO ARCS.
C
C                      I2 - LAST (LEFTMOST) BOUNDARY NODE
C                           WHICH IS VISIBLE FROM KK.
C
C                    IADJ - SET OF ADJACENCY LISTS OF NODES
C                           IN THE MESH.
C
C                    IEND - POINTERS TO THE ENDS OF
C                           ADJACENCY LISTS IN IADJ FOR
C                           EACH NODE IN THE MESH.
C
C   IADJ AND IEND MAY BE CREATED BY TRMESH AND MUST CONTAIN
C THE VERTICES I1 AND I2.  I1 AND I2 MAY BE DETERMINED BY
C TRFIND.
C
C KK, I1, AND I2 ARE NOT ALTERED BY THIS ROUTINE.
C
C OUTPUT PARAMETERS - IADJ,IEND - UPDATED WITH THE ADDITION
C                                 OF NODE KK AS THE LAST
C                                 ENTRY.  NODE KK WILL BE
C                                 CONNECTED TO I1, I2, AND
C                                 ALL BOUNDARY NODES BETWEEN
C                                 THEM.  NO OPTIMIZATION OF
C                                 THE MESH IS PERFORMED.
C
C MODULE REFERENCED BY BDYADD - SHIFTD
C
C INTRINSIC FUNCTIONS CALLED BY BDYADD - MIN0, MAX0
C
C***********************************************************
C
      INTEGER K, KM1, NRIGHT, NLEFT, NF, NL, N1, N2, I,
     .        IMIN, IMAX, KEND, NEXT, INDX
C
C LOCAL PARAMETERS -
C
C K =            LOCAL COPY OF KK
C KM1 =          K - 1
C NRIGHT,NLEFT = LOCAL COPIES OF I1, I2
C NF,NL =        INDICES OF IADJ BOUNDING THE PORTION OF THE
C                  ARRAY TO BE SHIFTED
C N1 =           IADJ INDEX OF THE FIRST NEIGHBOR OF NLEFT
C N2 =           IADJ INDEX OF THE LAST NEIGHBOR OF NRIGHT
C I =            DO-LOOP INDEX
C IMIN,IMAX =    BOUNDS ON DO-LOOP INDEX -- FIRST AND LAST
C                  ELEMENTS OF IEND TO BE INCREMENTED
C KEND =         POINTER TO THE LAST NEIGHBOR OF K IN IADJ
C NEXT =         NEXT BOUNDARY NODE TO BE CONNECTED TO KK
C INDX =         INDEX FOR IADJ
C
      K = KK
      KM1 = K - 1
      NRIGHT = I1
      NLEFT = I2
C
C INITIALIZE VARIABLES
C
      NL = IEND(KM1)
      N1 = 1
      IF (NLEFT .NE. 1) N1 = IEND(NLEFT-1) + 1
      N2 = IEND(NRIGHT)
      NF = MAX0(N1,N2)
C
C INSERT K AS A NEIGHBOR OF MAX(NRIGHT,NLEFT)
C
      CALL SHIFTD(NF,NL,2, IADJ)
      IADJ(NF+1) = K
      IMIN = MAX0(NRIGHT,NLEFT)
      DO 1 I = IMIN,KM1
    1   IEND(I) = IEND(I) + 2
C
C INITIALIZE KEND AND INSERT K AS A NEIGHBOR OF
C   MIN(NRIGHT,NLEFT)
C
      KEND = NL + 3
      NL = NF - 1
      NF = MIN0(N1,N2)
      CALL SHIFTD(NF,NL,1, IADJ)
      IADJ(NF) = K
      IMAX = IMIN - 1
      IMIN = MIN0(NRIGHT,NLEFT)
      DO 2 I = IMIN,IMAX
    2   IEND(I) = IEND(I) + 1
C
C INSERT NRIGHT AS THE FIRST NEIGHBOR OF K
C
      IADJ(KEND) = NRIGHT
C
C INITIALIZE INDX FOR LOOP ON BOUNDARY NODES BETWEEN NRIGHT
C   AND NLEFT
C
      INDX = IEND(NRIGHT) - 2
    3 NEXT = IADJ(INDX)
      IF (NEXT .EQ. NLEFT) GO TO 4
C
C CONNECT NEXT AND K
C
      KEND = KEND + 1
      IADJ(KEND) = NEXT
      INDX = IEND(NEXT)
      IADJ(INDX) = K
      INDX = INDX - 1
      GO TO 3
C
C INSERT NLEFT AND 0 AS THE LAST NEIGHBORS OF K
C
    4 IADJ(KEND+1) = NLEFT
      KEND = KEND + 2
      IADJ(KEND) = 0
      IEND(K) = KEND
      RETURN
      END
      SUBROUTINE INTADD (KK,I1,I2,I3, IADJ,IEND )
      INTEGER KK, I1, I2, I3, IADJ(*), IEND(KK)
C
C***********************************************************
C
C                                               ROBERT RENKA
C                                       OAK RIDGE NATL. LAB.
C                                             (615) 576-5139
C
C   THIS ROUTINE ADDS AN INTERIOR NODE TO A TRIANGULATION
C OF A SET OF KK-1 POINTS ON THE UNIT SPHERE.  IADJ AND IEND
C ARE UPDATED WITH THE INSERTION OF NODE KK IN THE TRIANGLE
C WHOSE VERTICES ARE I1, I2, AND I3.
C
C INPUT PARAMETERS -        KK - INDEX OF NODE TO BE
C                                INSERTED.  KK .GE. 4.
C
C                     I1,I2,I3 - INDICES OF THE VERTICES OF
C                                A TRIANGLE CONTAINING NODE
C                                KK - IN COUNTERCLOCKWISE
C                                ORDER.
C
C                         IADJ - SET OF ADJACENCY LISTS
C                                OF NODES IN THE MESH.
C
C                         IEND - POINTERS TO THE ENDS OF
C                                ADJACENCY LISTS IN IADJ FOR
C                                EACH NODE IN THE MESH.
C
C   IADJ AND IEND MAY BE CREATED BY TRMESH AND MUST CONTAIN
C THE VERTICES I1, I2, AND I3.  I1,I2,I3 MAY BE DETERMINED
C BY TRFIND.
C
C KK, I1, I2, AND I3 ARE NOT ALTERED BY THIS ROUTINE.
C
C OUTPUT PARAMETERS - IADJ,IEND - UPDATED WITH THE ADDITION
C                                 OF NODE KK AS THE LAST
C                                 ENTRY.  NODE KK WILL BE
C                                 CONNECTED TO NODES I1, I2,
C                                 AND I3.  NO OPTIMIZATION
C                                 OF THE MESH IS PERFORMED.
C
C MODULE REFERENCED BY INTADD - SHIFTD
C
C INTRINSIC FUNCTION CALLED BY INTADD - MOD
C
C***********************************************************
C
      INTEGER K, KM1, N(3), NFT(3), IP1, IP2, IP3, INDX, NF,
     .        NL, N1, N2, IMIN, IMAX, I, ITEMP
C
C LOCAL PARAMETERS -
C
C K =           LOCAL COPY OF KK
C KM1 =         K - 1
C N =           VECTOR CONTAINING I1, I2, I3
C NFT =         POINTERS TO THE TOPS OF THE 3 SETS OF IADJ
C                 ELEMENTS TO BE SHIFTED DOWNWARD
C IP1,IP2,IP3 = PERMUTATION INDICES FOR N AND NFT
C INDX =        INDEX FOR IADJ AND N
C NF,NL =       INDICES OF FIRST AND LAST ENTRIES IN IADJ
C                 TO BE SHIFTED DOWN
C N1,N2 =       FIRST 2 VERTICES OF A NEW TRIANGLE --
C                 (N1,N2,KK)
C IMIN,IMAX =   BOUNDS ON DO-LOOP INDEX -- FIRST AND LAST
C                 ELEMENTS OF IEND TO BE INCREMENTED
C I =           DO-LOOP INDEX
C ITEMP =       TEMPORARY STORAGE LOCATION
C
      K = KK
C
C INITIALIZATION
C
      N(1) = I1
      N(2) = I2
      N(3) = I3
C
C SET UP NFT
C
      DO 2 I = 1,3
        N1 = N(I)
        INDX = MOD(I,3) + 1
        N2 = N(INDX)
        INDX = IEND(N1) + 1
C
C FIND THE INDEX OF N2 AS A NEIGHBOR OF N1
C
    1   INDX = INDX - 1
        IF (IADJ(INDX) .NE. N2) GO TO 1
    2   NFT(I) = INDX + 1
C
C ORDER THE VERTICES BY DECREASING MAGNITUDE -
C   N(IP(I+1)) PRECEDES N(IP(I)) IN IEND FOR
C   I = 1,2.
C
      IP1 = 1
      IP2 = 2
      IP3 = 3
      IF ( N(2) .LE. N(1) ) GO TO 3
      IP1 = 2
      IP2 = 1
    3 IF ( N(3) .LE. N(IP1) ) GO TO 4
      IP3 = IP1
      IP1 = 3
    4 IF ( N(IP3) .LE. N(IP2) )  GO TO 5
      ITEMP = IP2
      IP2 = IP3
      IP3 = ITEMP
C
C ADD NODE K TO THE ADJACENCY LISTS OF EACH VERTEX AND
C   UPDATE IEND.  FOR EACH VERTEX, A SET OF IADJ ELEMENTS
C   IS SHIFTED DOWNWARD AND K IS INSERTED.  SHIFTING STARTS
C   AT THE END OF THE ARRAY.
C
    5 KM1 = K - 1
      NL = IEND(KM1)
      NF = NFT(IP1)
      IF (NF .LE. NL) CALL SHIFTD(NF,NL,3, IADJ)
      IADJ(NF+2) = K
      IMIN = N(IP1)
      IMAX = KM1
      DO 6 I = IMIN,IMAX
    6   IEND(I) = IEND(I) + 3
C
      NL = NF - 1
      NF = NFT(IP2)
      CALL SHIFTD(NF,NL,2, IADJ)
      IADJ(NF+1) = K
      IMAX = IMIN - 1
      IMIN = N(IP2)
      DO 7 I = IMIN,IMAX
    7   IEND(I) = IEND(I) + 2
C
      NL = NF - 1
      NF = NFT(IP3)
      CALL SHIFTD(NF,NL,1, IADJ)
      IADJ(NF) = K
      IMAX = IMIN - 1
      IMIN = N(IP3)
      DO 8 I = IMIN,IMAX
    8   IEND(I) = IEND(I) + 1
C
C ADD NODE K TO IEND AND ITS NEIGHBORS TO IADJ
C
      INDX = IEND(KM1)
      IEND(K) = INDX + 3
      DO 9 I = 1,3
        INDX = INDX + 1
    9   IADJ(INDX) = N(I)
      RETURN
      END
      SUBROUTINE SWAP (NIN1,NIN2,NOUT1,NOUT2, IADJ,IEND )
      INTEGER NIN1, NIN2, NOUT1, NOUT2, IADJ(*), IEND(*)
      EXTERNAL INDEX
C
C***********************************************************
C
C                                               ROBERT RENKA
C                                       OAK RIDGE NATL. LAB.
C                                             (615) 576-5139
C
C   THIS SUBROUTINE SWAPS THE DIAGONALS IN A CONVEX QUADRI-
C LATERAL.
C
C INPUT PARAMETERS -  NIN1,NIN2,NOUT1,NOUT2 - NODAL INDICES
C                            OF A PAIR OF ADJACENT TRIANGLES
C                            WHICH FORM A CONVEX QUADRILAT-
C                            ERAL.  NOUT1 AND NOUT2 ARE CON-
C                            NECTED BY AN ARC WHICH IS TO BE
C                            REPLACED BY THE ARC NIN1-NIN2.
C                            (NIN1,NOUT1,NOUT2) MUST BE TRI-
C                            ANGLE VERTICES IN COUNTERCLOCK-
C                            WISE ORDER.
C
C THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
C
C                IADJ,IEND - TRIANGULATION DATA STRUCTURE
C                            (SEE SUBROUTINE TRMESH).
C
C OUTPUT PARAMETERS - IADJ,IEND - UPDATED WITH THE ARC
C                                 REPLACEMENT.
C
C MODULES REFERENCED BY SWAP - INDEX, SHIFTD
C
C***********************************************************
C
      INTEGER IN(2), IO(2), IP1, IP2, J, K, NF, NL, I,
     .        IMIN, IMAX
C
C LOCAL PARAMETERS -
C
C IN =        NIN1 AND NIN2 ORDERED BY INCREASING MAGNITUDE
C               (THE NEIGHBORS OF IN(1) PRECEDE THOSE OF
C               IN(2) IN IADJ)
C IO =        NOUT1 AND NOUT2 IN INCREASING ORDER
C IP1,IP2 =   PERMUTATION OF (1,2) SUCH THAT IO(IP1)
C               PRECEDES IO(IP2) AS A NEIGHBOR OF IN(1)
C J,K =       PERMUTATION OF (1,2) USED AS INDICES OF IN
C               AND IO
C NF,NL =     IADJ INDICES BOUNDARY A PORTION OF THE ARRAY
C               TO BE SHIFTED
C I =         IEND INDEX
C IMIN,IMAX = BOUNDS ON THE PORTION OF IEND TO BE INCRE-
C               MENTED OR DECREMENTED
C
      IN(1) = NIN1
      IN(2) = NIN2
      IO(1) = NOUT1
      IO(2) = NOUT2
      IP1 = 1
C
C ORDER THE INDICES SO THAT IN(1) .LT. IN(2) AND IO(1) .LT.
C   IO(2), AND CHOOSE IP1 AND IP2 SUCH THAT (IN(1),IO(IP1),
C   IO(IP2)) FORMS A TRIANGLE.
C
      IF (IN(1) .LT. IN(2)) GO TO 1
      IN(1) = IN(2)
      IN(2) = NIN1
      IP1 = 2
    1 IF (IO(1) .LT. IO(2)) GO TO 2
      IO(1) = IO(2)
      IO(2) = NOUT1
      IP1 = 3 - IP1
    2 IP2 = 3 - IP1
      IF (IO(2) .LT. IN(1)) GO TO 8
      IF (IN(2) .LT. IO(1)) GO TO 12
C
C IN(1) AND IO(1) PRECEDE IN(2) AND IO(2).  FOR (J,K) =
C   (1,2) AND (2,1), DELETE IO(K) AS A NEIGHBOR OF IO(J)
C   BY SHIFTING A PORTION OF IADJ EITHER UP OR DOWN AND
C   AND INSERT IN(K) AS A NEIGHBOR OF IN(J).
C
      DO 7 J = 1,2
        K = 3 - J
        IF (IN(J) .GT. IO(J)) GO TO 4
C
C   THE NEIGHBORS OF IN(J) PRECEDE THOSE OF IO(J) -- SHIFT
C     DOWN BY 1
C
        NF = 1 + INDEX(IN(J),IO(IP1),IADJ,IEND)
        NL = -1 + INDEX(IO(J),IO(K),IADJ,IEND)
        IF (NF .LE. NL) CALL SHIFTD(NF,NL,1, IADJ )
        IADJ(NF) = IN(K)
        IMIN = IN(J)
        IMAX = IO(J)-1
        DO 3 I = IMIN,IMAX
    3     IEND(I) = IEND(I) + 1
        GO TO 6
C
C   THE NEIGHBORS OF IO(J) PRECEDE THOSE OF IN(J) -- SHIFT
C     UP BY 1
C
    4   NF = 1 + INDEX(IO(J),IO(K),IADJ,IEND)
        NL = -1 + INDEX(IN(J),IO(IP2),IADJ,IEND)
        IF (NF .LE. NL) CALL SHIFTD(NF,NL,-1, IADJ )
        IADJ(NL) = IN(K)
        IMIN = IO(J)
        IMAX = IN(J) - 1
        DO 5 I = IMIN,IMAX
    5     IEND(I) = IEND(I) - 1
C
C   REVERSE (IP1,IP2) FOR (J,K) = (2,1)
C
    6   IP1 = IP2
        IP2 = 3 - IP1
    7   CONTINUE
      RETURN
C
C THE VERTICES ARE ORDERED (IO(1),IO(2),IN(1),IN(2)).
C   DELETE IO(2) BY SHIFTING UP BY 1
C
    8 NF = 1 + INDEX(IO(1),IO(2),IADJ,IEND)
      NL = -1 + INDEX(IO(2),IO(1),IADJ,IEND)
      IF (NF .LE. NL) CALL SHIFTD(NF,NL,-1, IADJ )
      IMIN = IO(1)
      IMAX = IO(2)-1
      DO 9 I = IMIN,IMAX
    9   IEND(I) = IEND(I) - 1
C
C   DELETE IO(1) BY SHIFTING UP BY 2 AND INSERT IN(2)
C
      NF = NL + 2
      NL = -1 + INDEX(IN(1),IO(IP2),IADJ,IEND)
      IF (NF .LE. NL) CALL SHIFTD(NF,NL,-2, IADJ )
      IADJ(NL-1) = IN(2)
      IMIN = IO(2)
      IMAX = IN(1)-1
      DO 10 I = IMIN,IMAX
   10   IEND(I) = IEND(I) - 2
C
C   SHIFT UP BY 1 AND INSERT IN(1)
C
      NF = NL + 1
      NL = -1 + INDEX(IN(2),IO(IP1),IADJ,IEND)
      CALL SHIFTD(NF,NL,-1, IADJ )
      IADJ(NL) = IN(1)
      IMIN = IN(1)
      IMAX = IN(2)-1
      DO 11 I = IMIN,IMAX
   11   IEND(I) = IEND(I) - 1
      RETURN
C
C THE VERTICES ARE ORDERED (IN(1),IN(2),IO(1),IO(2)).
C   DELETE IO(1) BY SHIFTING DOWN BY 1
C
   12 NF = 1 + INDEX(IO(1),IO(2),IADJ,IEND)
      NL = -1 + INDEX(IO(2),IO(1),IADJ,IEND)
      IF (NF .LE. NL) CALL SHIFTD(NF,NL,1, IADJ )
      IMIN = IO(1)
      IMAX = IO(2) - 1
      DO 13 I = IMIN,IMAX
   13   IEND(I) = IEND(I) + 1
C
C   DELETE IO(2) BY SHIFTING DOWN BY 2 AND INSERT IN(1)
C
      NL = NF - 2
      NF = 1 + INDEX(IN(2),IO(IP2),IADJ,IEND)
      IF (NF .LE. NL) CALL SHIFTD(NF,NL,2, IADJ )
      IADJ(NF+1) = IN(1)
      IMIN = IN(2)
      IMAX = IO(1) - 1
      DO 14 I = IMIN,IMAX
   14   IEND(I) = IEND(I) + 2
C
C   SHIFT DOWN BY 1 AND INSERT IN(2)
C
      NL = NF - 1
      NF = 1 + INDEX(IN(1),IO(IP1),IADJ,IEND)
      CALL SHIFTD(NF,NL,1, IADJ )
      IADJ(NF) = IN(2)
      IMIN = IN(1)
      IMAX = IN(2) - 1
      DO 15 I = IMIN,IMAX
   15   IEND(I) = IEND(I) + 1
      RETURN
      END
      SUBROUTINE TRFIND (NST,PX,PY,X,Y,IADJ,IEND, I1,I2,I3)
      INTEGER NST, IADJ(*), IEND(*), I1, I2, I3
      REAL    PX, PY, X(*), Y(*)
C
C***********************************************************
C
C                                               ROBERT RENKA
C                                       OAK RIDGE NATL. LAB.
C                                             (615) 576-5139
C
C   THIS ROUTINE LOCATES A POINT P IN A THIESSEN TRIANGU-
C LATION, RETURNING THE VERTEX INDICES OF A TRIANGLE WHICH
C CONTAINS P.  TRFIND IS PART OF AN INTERPOLATION PACKAGE
C WHICH PROVIDES SUBROUTINES FOR CREATING THE MESH.
C
C INPUT PARAMETERS -    NST - INDEX OF NODE AT WHICH TRFIND
C                             BEGINS SEARCH.  SEARCH TIME
C                             DEPENDS ON THE PROXIMITY OF
C                             NST TO P.
C
C                     PX,PY - X AND Y-COORDINATES OF THE
C                             POINT TO BE LOCATED.
C
C                       X,Y - VECTORS OF COORDINATES OF
C                             NODES IN THE MESH.  (X(I),Y(I))
C                             DEFINES NODE I FOR I = 1,...,N
C                             WHERE N .GE. 3.
C
C                      IADJ - SET OF ADJACENCY LISTS OF
C                             NODES IN THE MESH.
C
C                      IEND - POINTERS TO THE ENDS OF
C                             ADJACENCY LISTS IN IADJ FOR
C                             EACH NODE IN THE MESH.
C
C IADJ AND IEND MAY BE CREATED BY TRMESH.
C
C INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
C
C OUTPUT PARAMETERS - I1,I2,I3 - VERTEX INDICES IN COUNTER-
C                                CLOCKWISE ORDER - VERTICES
C                                OF A TRIANGLE CONTAINING P
C                                IF P IS AN INTERIOR NODE.
C                                IF P IS OUTSIDE OF THE
C                                BOUNDARY OF THE MESH, I1
C                                AND I2 ARE THE FIRST (RIGHT
C                                -MOST) AND LAST (LEFTMOST)
C                                NODES WHICH ARE VISIBLE
C                                FROM P, AND I3 = 0.  IF P
C                                AND ALL OF THE NODES LIE ON
C                                A SINGLE LINE THEN I1 = I2
C                                = I3 = 0.
C
C MODULES REFERENCED BY TRFIND - NONE
C
C INTRINSIC FUNCTION CALLED BY TRFIND - MAX0
C
C***********************************************************
C
      INTEGER N0, N1, N2, N3, N4, INDX, IND, NF,
     .        NL, NEXT
      REAL    XP, YP
      LOGICAL LEFT
C
C LOCAL PARAMETERS -
C
C XP,YP =     LOCAL VARIABLES CONTAINING PX AND PY
C N0,N1,N2 =  NODES IN COUNTERCLOCKWISE ORDER DEFINING A
C               CONE (WITH VERTEX N0) CONTAINING P
C N3,N4 =     NODES OPPOSITE N1-N2 AND N2-N1, RESPECTIVELY
C INDX,IND =  INDICES FOR IADJ
C NF,NL =     FIRST AND LAST NEIGHBORS OF N0 IN IADJ, OR
C               FIRST (RIGHTMOST) AND LAST (LEFTMOST) NODES
C               VISIBLE FROM P WHEN P IS OUTSIDE THE
C               BOUNDARY
C NEXT =      CANDIDATE FOR I1 OR I2 WHEN P IS OUTSIDE OF
C               THE BOUNDARY
C LEFT =      STATEMENT FUNCTION WHICH COMPUTES THE SIGN OF
C               A CROSS PRODUCT (Z-COMPONENT).  LEFT(X1,...,
C               Y0) = .TRUE. IFF (X0,Y0) IS ON OR TO THE
C               LEFT OF THE VECTOR FROM (X1,Y1) TO (X2,Y2).
C
      LEFT(X1,Y1,X2,Y2,X0,Y0) = (X2-X1)*(Y0-Y1) .GE.
     .                          (X0-X1)*(Y2-Y1)
      XP = PX
      YP = PY
C
C INITIALIZE VARIABLES AND FIND A CONE CONTAINING P
C
      N0 = MAX0(NST,1)
    1 INDX = IEND(N0)
      NL = IADJ(INDX)
      INDX = 1
      IF (N0 .NE. 1) INDX = IEND(N0-1) + 1
      NF = IADJ(INDX)
      N1 = NF
      IF (NL .NE. 0) GO TO 3
C
C N0 IS A BOUNDARY NODE.  SET NL TO THE LAST NONZERO
C   NEIGHBOR OF N0.
C
      IND = IEND(N0) - 1
      NL = IADJ(IND)
      IF ( LEFT(X(N0),Y(N0),X(NF),Y(NF),XP,YP) ) GO TO 2
C
C P IS OUTSIDE THE BOUNDARY
C
      NL = N0
      GO TO 16
    2 IF ( LEFT(X(NL),Y(NL),X(N0),Y(N0),XP,YP) ) GO TO 4
C
C P IS OUTSIDE THE BOUNDARY AND N0 IS THE RIGHTMOST
C   VISIBLE BOUNDARY NODE
C
      I1 = N0
      GO TO 18
C
C N0 IS AN INTERIOR NODE.  FIND N1.
C
    3 IF ( LEFT(X(N0),Y(N0),X(N1),Y(N1),XP,YP) ) GO TO 4
      INDX = INDX + 1
      N1 = IADJ(INDX)
      IF (N1 .EQ. NL) GO TO 7
      GO TO 3
C
C P IS TO THE LEFT OF ARC N0-N1.  INITIALIZE N2 TO THE NEXT
C   NEIGHBOR OF N0.
C
    4 INDX = INDX + 1
      N2 = IADJ(INDX)
      IF ( .NOT. LEFT(X(N0),Y(N0),X(N2),Y(N2),XP,YP) )
     .   GO TO 8
      N1 = N2
      IF (N1 .NE. NL) GO TO 4
      IF ( .NOT. LEFT(X(N0),Y(N0),X(NF),Y(NF),XP,YP) )
     .   GO TO 7
      IF (XP .EQ. X(N0) .AND. YP .EQ. Y(N0)) GO TO 6
C
C P IS LEFT OF OR ON ARCS N0-NB FOR ALL NEIGHBORS NB
C   OF N0.
C ALL POINTS ARE COLLINEAR IFF P IS LEFT OF NB-N0 FOR
C   ALL NEIGHBORS NB OF N0.  SEARCH THE NEIGHBORS OF N0
C   IN REVERSE ORDER.  NOTE -- N1 = NL AND INDX POINTS TO
C   NL.
C
    5 IF ( .NOT. LEFT(X(N1),Y(N1),X(N0),Y(N0),XP,YP) )
     .   GO TO 6
      IF (N1 .EQ. NF) GO TO 20
      INDX = INDX - 1
      N1 = IADJ(INDX)
      GO TO 5
C
C P IS TO THE RIGHT OF N1-N0, OR P=N0.  SET N0 TO N1 AND
C   START OVER.
C
    6 N0 = N1
      GO TO 1
C
C P IS BETWEEN ARCS N0-N1 AND N0-NF
C
    7 N2 = NF
C
C P IS CONTAINED IN A CONE DEFINED BY LINE SEGMENTS N0-N1
C   AND N0-N2 WHERE N1 IS ADJACENT TO N2
C
    8 N3 = N0
    9 IF ( LEFT(X(N1),Y(N1),X(N2),Y(N2),XP,YP) ) GO TO 13
C
C SET N4 TO THE FIRST NEIGHBOR OF N2 FOLLOWING N1
C
      INDX = IEND(N2)
      IF (IADJ(INDX) .NE. N1) GO TO 10
C
C N1 IS THE LAST NEIGHBOR OF N2.
C SET N4 TO THE FIRST NEIGHBOR.
C
      INDX = 1
      IF (N2 .NE. 1) INDX = IEND(N2-1) + 1
      N4 = IADJ(INDX)
      GO TO 11
C
C N1 IS NOT THE LAST NEIGHBOR OF N2
C
   10 INDX = INDX-1
      IF (IADJ(INDX) .NE. N1) GO TO 10
      N4 = IADJ(INDX+1)
      IF (N4 .NE. 0) GO TO 11
C
C P IS OUTSIDE THE BOUNDARY
C
      NF = N2
      NL = N1
      GO TO 16
C
C DEFINE A NEW ARC N1-N2 WHICH INTERSECTS THE LINE
C   SEGMENT N0-P
C
   11 IF ( LEFT(X(N0),Y(N0),X(N4),Y(N4),XP,YP) ) GO TO 12
      N3 = N2
      N2 = N4
      GO TO 9
   12 N3 = N1
      N1 = N4
      GO TO 9
C
C P IS IN THE TRIANGLE (N1,N2,N3) AND NOT ON N2-N3.  IF
C   N3-N1 OR N1-N2 IS A BOUNDARY ARC CONTAINING P, TREAT P
C   AS EXTERIOR.
C
   13 INDX = IEND(N1)
      IF (IADJ(INDX) .NE. 0) GO TO 15
C
C N1 IS A BOUNDARY NODE.  N3-N1 IS A BOUNDARY ARC IFF N3
C   IS THE LAST NONZERO NEIGHBOR OF N1.
C
      IF (N3 .NE. IADJ(INDX-1)) GO TO 14
C
C N3-N1 IS A BOUNDARY ARC
C
      IF ( .NOT. LEFT(X(N1),Y(N1),X(N3),Y(N3),XP,YP) )
     .   GO TO 14
C
C P LIES ON N1-N3
C
      I1 = N1
      I2 = N3
      I3 = 0
      RETURN
C
C N3-N1 IS NOT A BOUNDARY ARC CONTAINING P.  N1-N2 IS A
C   BOUNDARY ARC IFF N2 IS THE FIRST NEIGHBOR OF N1.
C
   14 INDX = 1
      IF (N1 .NE. 1) INDX = IEND(N1-1) + 1
      IF (N2 .NE. IADJ(INDX)) GO TO 15
C
C N1-N2 IS A BOUNDARY ARC
C
      IF ( .NOT. LEFT(X(N2),Y(N2),X(N1),Y(N1),XP,YP) )
     .   GO TO 15
C
C P LIES ON N1-N2
C
      I1 = N2
      I2 = N1
      I3 = 0
      RETURN
C
C P DOES NOT LIE ON A BOUNDARY ARC.
C
   15 I1 = N1
      I2 = N2
      I3 = N3
      RETURN
C
C NF AND NL ARE ADJACENT BOUNDARY NODES WHICH ARE VISIBLE
C   FROM P.  FIND THE FIRST VISIBLE BOUNDARY NODE.
C SET NEXT TO THE FIRST NEIGHBOR OF NF.
C
   16 INDX = 1
      IF (NF .NE. 1) INDX = IEND(NF-1) + 1
      NEXT = IADJ(INDX)
      IF ( LEFT(X(NF),Y(NF),X(NEXT),Y(NEXT),XP,YP) )
     .   GO TO 17
      NF = NEXT
      GO TO 16
C
C NF IS THE FIRST (RIGHTMOST) VISIBLE BOUNDARY NODE
C
   17 I1 = NF
C
C FIND THE LAST VISIBLE BOUNDARY NODE.  NL IS THE FIRST
C   CANDIDATE FOR I2.
C SET NEXT TO THE LAST NEIGHBOR OF NL.
C
   18 INDX = IEND(NL) - 1
      NEXT = IADJ(INDX)
      IF ( LEFT(X(NEXT),Y(NEXT),X(NL),Y(NL),XP,YP) )
     .   GO TO 19
      NL = NEXT
      GO TO 18
C
C NL IS THE LAST (LEFTMOST) VISIBLE BOUNDARY NODE
C
   19 I2 = NL
      I3 = 0
      RETURN
C
C ALL POINTS ARE COLLINEAR
C
   20 I1 = 0
      I2 = 0
      I3 = 0
      RETURN
      END
      LOGICAL FUNCTION SWPTST (IN1,IN2,IO1,IO2,X,Y)
      INTEGER IN1, IN2, IO1, IO2
      REAL    X(*), Y(*)
C
C***********************************************************
C
C                                               ROBERT RENKA
C                                       OAK RIDGE NATL. LAB.
C                                             (615) 576-5139
C
C   THIS FUNCTION DECIDES WHETHER OR NOT TO REPLACE A
C DIAGONAL ARC IN A QUADRILATERAL WITH THE OTHER DIAGONAL.
C THE DETERMINATION IS BASED ON THE SIZES OF THE ANGLES
C CONTAINED IN THE 2 TRIANGLES DEFINED BY THE DIAGONAL.
C THE DIAGONAL IS CHOSEN TO MAXIMIZE THE SMALLEST OF THE
C SIX ANGLES OVER THE TWO PAIRS OF TRIANGLES.
C
C INPUT PARAMETERS -  IN1,IN2,IO1,IO2 - NODE INDICES OF THE
C                              FOUR POINTS DEFINING THE
C                              QUADRILATERAL.  IO1 AND IO2
C                              ARE CURRENTLY CONNECTED BY A
C                              DIAGONAL ARC.  THIS ARC
C                              SHOULD BE REPLACED BY AN ARC
C                              CONNECTING IN1, IN2 IF THE
C                              DECISION IS MADE TO SWAP.
C                              IN1,IO1,IO2 MUST BE IN
C                              COUNTERCLOCKWISE ORDER.
C
C                        X,Y - VECTORS OF NODAL COORDINATES.
C                              (X(I),Y(I)) ARE THE COORD-
C                              INATES OF NODE I FOR I = IN1,
C                              IN2, IO1, OR IO2.
C
C NONE OF THE INPUT PARAMETERS ARE ALTERED BY THIS ROUTINE.
C
C OUTPUT PARAMETER -  SWPTST - .TRUE. IFF THE ARC CONNECTING
C                              IO1 AND IO2 IS TO BE REPLACED
C
C MODULES REFERENCED BY SWPTST - NONE
C
C***********************************************************
C
      REAL DX11, DX12, DX22, DX21, DY11, DY12, DY22, DY21,
     .     SIN1, SIN2, COS1, COS2, SIN12
C
C LOCAL PARAMETERS -
C
C DX11,DY11 = X,Y COORDINATES OF THE VECTOR IN1-IO1
C DX12,DY12 = X,Y COORDINATES OF THE VECTOR IN1-IO2
C DX22,DY22 = X,Y COORDINATES OF THE VECTOR IN2-IO2
C DX21,DY21 = X,Y COORDINATES OF THE VECTOR IN2-IO1
C SIN1 =      CROSS PRODUCT OF THE VECTORS IN1-IO1 AND
C               IN1-IO2 -- PROPORTIONAL TO SIN(T1) WHERE T1
C               IS THE ANGLE AT IN1 FORMED BY THE VECTORS
C COS1 =      INNER PRODUCT OF THE VECTORS IN1-IO1 AND
C               IN1-IO2 -- PROPORTIONAL TO COS(T1)
C SIN2 =      CROSS PRODUCT OF THE VECTORS IN2-IO2 AND
C               IN2-IO1 -- PROPORTIONAL TO SIN(T2) WHERE T2
C               IS THE ANGLE AT IN2 FORMED BY THE VECTORS
C COS2 =      INNER PRODUCT OF THE VECTORS IN2-IO2 AND
C               IN2-IO1 -- PROPORTIONAL TO COS(T2)
C SIN12 =     SIN1*COS2 + COS1*SIN2 -- PROPORTIONAL TO
C               SIN(T1+T2)
C
      SWPTST = .FALSE.
C
C COMPUTE THE VECTORS CONTAINING THE ANGLES T1, T2
C
      DX11 = X(IO1) - X(IN1)
      DX12 = X(IO2) - X(IN1)
      DX22 = X(IO2) - X(IN2)
      DX21 = X(IO1) - X(IN2)
C
      DY11 = Y(IO1) - Y(IN1)
      DY12 = Y(IO2) - Y(IN1)
      DY22 = Y(IO2) - Y(IN2)
      DY21 = Y(IO1) - Y(IN2)
C
C COMPUTE INNER PRODUCTS
C
      COS1 = DX11*DX12 + DY11*DY12
      COS2 = DX22*DX21 + DY22*DY21
C
C THE DIAGONALS SHOULD BE SWAPPED IFF (T1+T2) .GT. 180
C   DEGREES.  THE FOLLOWING TWO TESTS INSURE NUMERICAL
C   STABILITY.
C
      IF (COS1 .GE. 0.  .AND.  COS2 .GE. 0.) RETURN
      IF (COS1 .LT. 0.  .AND.  COS2 .LT. 0.) GO TO 1
C
C COMPUTE VECTOR CROSS PRODUCTS
C
      SIN1 = DX11*DY12 - DX12*DY11
      SIN2 = DX22*DY21 - DX21*DY22
      SIN12 = SIN1*COS2 + COS1*SIN2
      IF (SIN12 .GE. 0.) RETURN
    1 SWPTST = .TRUE.
      RETURN
      END
      INTEGER FUNCTION INDEX (NVERTX,NABOR,IADJ,IEND)
      INTEGER NVERTX, NABOR, IADJ(*), IEND(*)
C
C***********************************************************
C
C                                               ROBERT RENKA
C                                       OAK RIDGE NATL. LAB.
C                                             (615) 576-5139
C
C   THIS FUNCTION RETURNS THE INDEX OF NABOR IN THE
C ADJACENCY LIST FOR NVERTX.
C
C INPUT PARAMETERS - NVERTX - NODE WHOSE ADJACENCY LIST IS
C                             TO BE SEARCHED.
C
C                     NABOR - NODE WHOSE INDEX IS TO BE
C                             RETURNED.  NABOR MUST BE
C                             CONNECTED TO NVERTX.
C
C                      IADJ - SET OF ADJACENCY LISTS.
C
C                      IEND - POINTERS TO THE ENDS OF
C                             ADJACENCY LISTS IN IADJ.
C
C INPUT PARAMETERS ARE NOT ALTERED BY THIS FUNCTION.
C
C OUTPUT PARAMETER -  INDEX - IADJ(INDEX) = NABOR.
C
C MODULES REFERENCED BY INDEX - NONE
C
C***********************************************************
C
      INTEGER NB, INDX
C
C LOCAL PARAMETERS -
C
C NB =   LOCAL COPY OF NABOR
C INDX = INDEX FOR IADJ
C
      NB = NABOR
C
C INITIALIZATION
C
      INDX = IEND(NVERTX) + 1
C
C SEARCH THE LIST OF NVERTX NEIGHBORS FOR NB
C
    1 INDX = INDX - 1
      IF (IADJ(INDX) .NE. NB) GO TO 1
C
      INDEX = INDX
      RETURN
      END
