C**** Incorporates changes given in Remark by Einarsson
C**** CACM 17(4) April 1974 p.225
      REAL FUNCTION WEW_A ( X, EN )
C
C  ITERATIVE SOLUTION OF X = W * EXP ( W ) WHERE X IS GIVEN.
C  (NOVEMBER 1970)
C  (REVISED - SEPTEMBER 1971)
C  VERSION A -- CDC 6600 MACHINE ACCURACY.
C
C  INPUT PARAMETER:
C    X  ARGUMENT OF W(X)
C
C  OUTPUT PARAMETERS:
C    WEW  THE DESIRED SOLUTION.
C    EN   THE LAST RELATIVE CORRECTION TO W(X).
C
C  SET CONSTANTS...
C     .. Scalar Arguments ..
      REAL EN,X
C     ..
C     .. Local Scalars ..
C**** Next 2 statement replaced in remark by Einarsson
C     REAL C1,C2,C3,C4,FLOGX,TEMP,TEMP2,WN,Y,ZN
C     INTEGER NEWE
      REAL FLOGX,TEMP,TEMP2,WN,Y,ZN
C**** End of replaced statements
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC LOG
C**** Next 10 statements deleted in remark by Einarsson
C     ..
C     .. Save statement ..
C      SAVE C1,C2,C3,C4,NEWE
C      DATA NEWE / 1 /
C      IF ( NEWE ) 10, 20, 10
C10    NEWE = 0
C      C1 = 4. / 3.
C      C2 = 7. / 3.
C      C3 = 5. / 6.
C      C4 = 2. / 3.
C**** End of replaced statements
C
C  COMPUTE INITIAL GUESS...
20    FLOGX = LOG ( X )
      IF ( X - 6.46 ) 30, 30, 40
C**** Next statement replaced in remark by Einarsson
C30    WN = X * ( 1. + C1 * X ) / ( 1. + X * ( C2 + C3 * X ) )
30    WN = X*(3.0 + 4.0*X)/(3.0 + X*(7.0 + 2.5*X))
C**** End of replaced statements
      ZN = FLOGX - WN - LOG ( WN )
      GO TO 50
40    WN = FLOGX
      ZN = -LOG ( WN )
50    CONTINUE
C
C  ITERATION ONE...
      TEMP = 1. + WN
C**** Next statement replaced in remark by Einarsson
C     Y = 2. * TEMP * ( TEMP + C4 * ZN ) - ZN
      Y = 2. * TEMP * ( TEMP + ZN/1.5 ) - ZN
C**** End of replaced statements
      WN = WN * ( 1. + ZN * Y / ( TEMP * ( Y - ZN ) ) )
C
C  ITERATION TWO...
      ZN = FLOGX - WN - LOG ( WN )
      TEMP = 1. + WN
**** Next statement replaced in remark by Einarsson
C     TEMP2 = TEMP + C4 * ZN
      TEMP2 = TEMP + ZN/1.5
      EN = ZN * TEMP2 / ( TEMP * TEMP2 - .5 * ZN )
      WN = WN * ( 1. + EN )
C
C  RETURN...
      WEW_A = WN
      RETURN
      END
      REAL FUNCTION WEW_B ( X, EN )
C
C  ITERATIVE SOLUTION OF X = W * EXP ( W ) WHERE X IS GIVEN.
C  (NOVEMBER 1970)
C  (REVISED - SEPTEMBER 1971)
C  VERSION B -- MAXIMUM RELATIVE ERROR 3.E-7.
C
C  INPUT PARAMETER:
C    X  ARGUMENT OF W(X)
C
C  OUTPUT PARAMETERS:
C    WEW  THE DESIRED SOLUTION.
C    EN   THE LAST RELATIVE CORRECTION TO W(X).
C
C  SET CONSTANTS...
C     .. Scalar Arguments ..
      REAL EN,X,F
C     ..
C     .. Local Scalars ..
C**** Next 2 statement replaced in remark by Einarsson
C     REAL C1,C2,C3,C4,FLOGX,TEMP,WN,Y,ZN
C     INTEGER NEWE
      REAL FLOGX,TEMP,WN,Y,ZN
C**** End of replaced statements
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC LOG
C**** Next 10 statements deleted in remark by Einarsson
C     ..
C     .. Save statement ..
C      SAVE C1,C2,C3,C4,NEWE
C      DATA NEWE / 1 /
C      IF ( NEWE ) 10, 20, 10
C10    NEWE = 0
C      C1 = 4. / 3.
C      C2 = 7. / 3.
C      C3 = 5. / 6.
C      C4 = 2. / 3.
C**** End of replaced statements
C
C  COMPUTE INITIAL GUESS...
20    FLOGX = LOG ( X )
      IF ( X - .7385 ) 30, 30, 40
C**** Next statement replaced in remark by Einarsson
C30    WN = X * ( 1. + C1 * X ) / ( 1. + X * ( C2 + C3 * X ) )
30    WN = X*(3.0 + 4.0*X)/(3.0 + X*(7.0 + 2.5*X))
C**** End of replaced statements
      GO TO 50
40    WN = FLOGX - 24. * ((FLOGX+2.) * FLOGX-3.) /
     +     (( .7 * FLOGX + 58. ) * FLOGX + 127. )
50    CONTINUE
C
C  ITERATION ONE...
      ZN = FLOGX - WN - LOG ( WN )
      TEMP = 1. + WN
C**** Next statement replaced in remark by Einarsson
C     Y = 2. * TEMP * ( TEMP + C4 * ZN ) - ZN
      Y = 2. * TEMP * ( TEMP + ZN/1.5 ) - ZN
C**** End of replaced statements
      EN = ZN * Y / ( TEMP * ( Y - ZN ) )
      WN = WN * ( 1. + EN )
C
C  RETURN...
      WEW_B = WN
      RETURN
      END
