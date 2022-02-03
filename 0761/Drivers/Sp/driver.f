      PROGRAM  TPSD3P
C
C Test Program for the SDBI3P/SDSF3P subroutine package
C
C Hiroshi Akima
C U.S. Department of Commerce, NTIA/ITS
C Version of 1995/05
C
C This program calls the SDBI3P and SDSF3P subroutines.
C
C This program requires no input data files.
C
C This program creates the WFSD3P file.  All elements of the
C DZI1, DZI2, DZI3, and DZI4 arrays in this file are expected
C to be zero. 
C
C Specification Statements
      INTEGER NYI
      INTEGER NXI
      INTEGER I, IDP, IER, IWK, IXI, IYI, L, MD, NDP, NDPO2, NEWPG
      REAL DZI1, DZI2, DZI3, DZI4, WK, XD, XI, YD, YI, YIP, ZD, ZI1
     +, ZI2, ZI3, ZI4, ZIE
      PARAMETER   (NEWPG=210000000)
      PARAMETER   (NDP=30,NXI=6,NYI=5)
      CHARACTER   NMPR*6,NMWF*6
      DIMENSION   XD(NDP),YD(NDP),ZD(NDP),XI(NXI),YI(NYI),
     1            ZIE(NXI,NYI),
     2            ZI1(NXI,NYI),DZI1(NXI,NYI),
     3            ZI2(NXI,NYI),DZI2(NXI,NYI),
     4            ZI3(NXI,NYI),DZI3(NXI,NYI),
     5            ZI4(NXI,NYI),DZI4(NXI,NYI),
     6            YIP(NXI),
     7            WK(30,17),IWK(30,39)
C Data statements 
      DATA  NMPR/'TPSD3P'/, NMWF/'res'/
      DATA  (XD(I),YD(I),ZD(I),I=1,NDP)/
     1  11.16,  1.24, 22.15,   0.00,  0.00, 58.20,
     2  24.20, 16.23,  2.83,   9.66, 20.00,  4.73,
     3  19.85, 10.72,  7.97,   5.22, 14.66, 40.36,
     4  10.35,  4.11, 22.33,  11.77, 10.47, 13.62,
     5  19.72,  1.39, 16.83,  15.10, 17.19, 12.57,
     6   0.00, 20.00, 34.60,  25.00,  3.87,  8.74,
     7  20.87, 20.00,  5.74,  25.00,  0.00, 12.00,
     8  19.99,  4.62, 14.72,  14.59,  8.71, 14.81,
     9  10.28, 15.16, 21.59,  15.20,  0.00, 21.60,
     O   4.51, 20.00, 15.61,   5.23, 10.72, 26.50,
     1   0.00,  4.48, 61.77,   2.14, 15.03, 53.10,
     2  16.70, 19.65,  6.31,   0.51,  8.37, 49.43,
     3   6.08,  4.58, 35.74,  25.00, 20.00,  0.60,
     4  25.00, 11.87,  4.40,  21.67, 14.36,  5.52,
     5  14.90,  3.12, 21.70,   3.31,  0.13, 44.08/
      DATA  XI/
     1    0.00,  5.00, 10.00, 15.00, 20.00, 25.00/
      DATA  YI/ 
     1    0.00,  5.00, 10.00, 15.00, 20.00/ 
      DATA  ((ZIE(IXI,IYI),IYI=1,NYI),IXI=1,NXI)/ 
     1   58.20, 60.75, 49.33, 71.21, 34.60,
     2   37.28, 39.32, 27.61, 41.36, 14.83,
     3   25.29, 21.91, 16.36, 22.34,  4.12,
     4   21.69, 20.17, 12.98, 12.24,  3.73,
     5   17.56, 14.44,  8.23,  7.52,  6.32,
     6   12.00,  8.00,  5.27,  2.97,  0.60/
C Opens the output file.
      OPEN  (6,FILE=NMWF)
C Writes the input data.
      WRITE (6,6010)  NMPR,NDP
      NDPO2=NDP/2 
      DO 13  L=1,NDPO2 
        IF (MOD(L,5).EQ.1)   WRITE (6,6011) 
        WRITE (6,6012)
     1    (IDP,XD(IDP),YD(IDP),ZD(IDP),IDP=L,NDP,NDPO2) 
   13 CONTINUE 
C Calculates and writes the output results.
C Part 1.  Program check for SDBI3P
C   Subpart 1.1.  Calculation of each ZI value at a point
      DO 22  IYI=1,NYI 
        DO 21  IXI=1,NXI 
          IF (IXI.EQ.1.AND.IYI.EQ.1)  THEN
            MD=1
          ELSE
            MD=3
          END IF 
          CALL SDBI3P(MD,NDP,XD,YD,ZD,1,XI(IXI),YI(IYI), 
     1                ZI1(IXI,IYI),IER, WK,IWK)
          IF (IER.GT.0)  STOP 
          DZI1(IXI,IYI)=ABS(ZI1(IXI,IYI)-ZIE(IXI,IYI))
   21   CONTINUE 
   22 CONTINUE 
      WRITE (6,6025)  NEWPG,NMPR
      WRITE (6,6026)
      WRITE (6,6027)  YI,YI 
      WRITE (6,6028)  (XI(IXI),(ZI1(IXI,IYI),IYI=1,NYI),
     1                (DZI1(IXI,IYI),IYI=1,NYI),IXI=1,NXI)
C   Subpart 1.2.  Calculation of ZI values for each YI value 
      DO 33  IYI=1,NYI 
        IF (IYI.EQ.1)  THEN
          MD=1
        ELSE
          MD=3
        END IF 
        DO 31  IXI=1,NXI 
          YIP(IXI)=YI(IYI) 
   31   CONTINUE 
        CALL SDBI3P(MD,NDP,XD,YD,ZD,NXI,XI,YIP,
     1              ZI2(1,IYI),IER, WK,IWK)
        IF (IER.GT.0)  STOP 
        DO 32  IXI=1,NXI 
          DZI2(IXI,IYI)=ABS(ZI2(IXI,IYI)-ZIE(IXI,IYI))
   32   CONTINUE 
   33 CONTINUE 
      WRITE (6,6036)
      WRITE (6,6027)  YI,YI 
      WRITE (6,6028)  (XI(IXI),(ZI2(IXI,IYI),IYI=1,NYI),
     1                (DZI2(IXI,IYI),IYI=1,NYI),IXI=1,NXI)
C Part 2.  Program check for SDSF3P
C   Subpart 2.1.  Calculation with MD=1
      CALL SDSF3P(1,NDP,XD,YD,ZD,NXI,XI,NYI,YI, ZI3,IER, WK,IWK)
      IF (IER.GT.0)  STOP 
      DO 42  IYI=1,NYI 
        DO 41  IXI=1,NXI 
          DZI3(IXI,IYI)=ABS(ZI3(IXI,IYI)-ZIE(IXI,IYI))
   41   CONTINUE 
   42 CONTINUE 
      WRITE (6,6045)  NEWPG,NMPR
      WRITE (6,6046)
      WRITE (6,6027)  YI,YI 
      WRITE (6,6028)  (XI(IXI),(ZI3(IXI,IYI),IYI=1,NYI),
     1                (DZI3(IXI,IYI),IYI=1,NYI),IXI=1,NXI)
C   Subpart 2.2.  Calculation with MD=2
      CALL SDSF3P(2,NDP,XD,YD,ZD,NXI,XI,NYI,YI, ZI4,IER, WK,IWK)
      IF (IER.GT.0)  STOP 
      DO 52  IYI=1,NYI 
        DO 51  IXI=1,NXI 
          DZI4(IXI,IYI)=ABS(ZI4(IXI,IYI)-ZIE(IXI,IYI))
   51   CONTINUE 
   52 CONTINUE 
      WRITE (6,6056)
      WRITE (6,6027)  YI,YI 
      WRITE (6,6028)  (XI(IXI),(ZI4(IXI,IYI),IYI=1,NYI),
     1                (DZI4(IXI,IYI),IYI=1,NYI),IXI=1,NXI)
      STOP 
C Format Statements
 6010 FORMAT (A6////'   Input data        NDP =',I3//
     1  2('      I      XD     YD     ZD    ')) 
 6011 FORMAT (1X) 
 6012 FORMAT (2(5X,I2,2X,3F7.2,3X)) 
 6025 FORMAT (A1,A6,14X, 
     1  'Part 1.  Program Check for SDBI3P')
 6026 FORMAT (1X////'Calculation by points'//
     1  21X,'ZI values',28X,'Differences'// 
     2  21X,'ZI1(XI,YI)',27X,'DZI1(XI,YI)')
 6027 FORMAT ('   XI    YI=',34X,'YI='/ 
     1  5X,2(2X,5F7.2)/)
 6028 FORMAT (1X/F5.2,2X,5F7.2,2X,5F7.2)
 6036 FORMAT (1X////'Calculation by columns'// 
     1  21X,'ZI values',28X,'Differences'// 
     2  21X,'ZI2(XI,YI)',27X,'DZI2(XI,YI)')
 6045 FORMAT (A1,A6,14X, 
     1  'Part 2.  Program Check for SDSF3P')
 6046 FORMAT (1X////'Calculation with MD=1'//
     1  21X,'ZI values',28X,'Differences'// 
     2  21X,'ZI3(XI,YI)',27X,'DZI3(XI,YI)')
 6056 FORMAT (1X////'Calculation with MD=2'//
     1  21X,'ZI values',28X,'Differences'// 
     2  21X,'ZI4(XI,YI)',27X,'DZI4(XI,YI)')
      END
