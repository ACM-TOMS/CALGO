*DECK,SYSTXT
          IDENT  SYSTEXT
 SYSTEXT  TITLE  BLA LIBRARY SYSTEMS TEXT
          STEXT
 INFTN    SPACE  2,10
**        INFTN  NAME,NUM      PARAMETER CONVERSION
*                              FOR NON-RUN COMPILER
*         ENTRY:
*                NAME  =  ENTRY/EXIT NAME
*                NUM   =  NUMBER OF PARAMETERS
*
*                *F = 0   -  COMPASS
*                   = 1   -  RUN/MNF
*                   = 2   -  FTN
*
 INFTN    MACRO  NAME,NUM      SET UP FOR NON-RUN COMPILER
          LOCAL  M,RETURN
 .1       IFEQ   *F,2
 .2       IFNE   NUM,0
 I        DECMIC 0
 M        MIN    6,NUM
 .D       DUP    M
 I        DECMIC 'I'+1
          SB'I'  X1
          SA1    A1+1
 .D       ENDD
 .2       IFGE   NUM,7
 .D       DUP    NUM-6
 I        DECMIC 'I'+1
          BX6    X1
          SA6    NAME-NUM-2+'I'
 .D       ENDD
 .2       ENDIF
          SX6    A0
          SA6    SETA0
          EQ     RETURN
 SETA0    BSS    1
 RETURN   BSS    0
 .1       ENDIF
 INFTN    ENDM
 OUTFTN   SPACE  2,10
**        OUTFTN NAME          EXIT FOR NON-RUN COMPILER
*
*         ENTRY:
*                NAME  =  ENTRY/EXIT NAME
*
 OUTFTN   MACRO  NAME          EXIT FOR NON RUN COMPILER
 .1       IFEQ   *F,2
          SA1    SETA0
          SA0    X1
 .1       ENDIF
          EQ     NAME
 OUTFTN   ENDM
CALL      SPACE    2,10
**        CALL     SUBR,(ARG,...,ARG)  STANDARD RUN/FTN SUBROUTINE CALL.
*
*         ENTRY:   *F = 2, GENERATE FTN CALLING SEQUENCE.
*                     ' 2, GENERATE RUN CALLING SEQUENCE.
*                  SUBR = NAME OF THE SUBROUTINE.
*                  ARG = (OPTIONAL), ADDRESS OF ARGUMENT TO BE PASSED.
*                        <CONSTANT>, THE ADDRESS OF THE LITERAL VALUE
*                        <CONSTANT> IS PASSED.
*                      = (OMITTED), IF *F' 2 IRUN  THE CORRESPONDING B
*                        REGISTER OR ARGUMENT ADDRESS LOCATION IS NOT
*                        SET AND IS ASSUMED TO BE PRESET BY THE USER.
*                      = (OMITTED), IF *F = 2 IFTN  THE CORRESPONDING
*                        ARGUMENT ADDRESS LOCATION IS SET TO:
*                                  42/0LNULL,18/*+400000B.
* RUN     USES:    X - 1,6,7.      (IF MORE THAN 6 ARGUMENTS SPECIFIED)
*                  B - 1,..,N.     (IF N ARGUMENTS SPECIFIED AND N @ 7)
*                  B - ALL.        (IF MORE THAN 6 ARGUMENTS SPECIFIED)
*                  A - 1,6,7.      (IF MORE THAN 6 ARGUMENTS SPECIFIED)
* FTN     USES:    X - 1.
*                  X - 1,6,7.      (IF AN ARGUMENT EXPRESSION CONTAINS
*                                   A REGISTER)
*                  B - NONE.
*                  A - 0,1.
*                  A - 0,1,6,7.    (IF AN ARGUMENT EXPRESSION CONTAINS
*                                   A REGISTER)
*         CALLS:   SUBR.
*         NOTE:    TRACE BACK INFORMATION IS NOT DEFINED UNLESS THE
*                  MICRO 'ENTRY' IS DEFINED. THE MICRO 'ENTRY' IS
*                  DEFINED WITHIN THE BEGIN MACRO. IF NO ARGUMENTS ARE
*                  SPECIFIED THEN THE COMMA AND PARENTHESES MAY BE
*                  OMITTED. BOTH RUN AND FTN STYLE SUBROUTINES DO NOT
*                  PRESERVE REGISTER CONTENTS, EXCEPT FTN SYTLE
*                  SUBROUTINES PRESERVE REGISTER A0. IN THE FTN CALLING
*                  SEQUENCE THE CONTENTS OF THE CALLER/S REGISTER A1 IS
*                  PRESERVED BY ENTERING IT INTO REGISTER A0 BEFORE THE
*                  SUBROUTINE IS ENTERED AND RESETTING REGISTER A1 TO
*                  THE PRESERVED REGISTER A0 ON RETURN FROM THE
*                  SUBROUTINE.
*
CALL      MACRO    SUBR,ARGS
.1        IFEQ     *F,2
          FTN=1    SUBR,(ARGS)
.1        ELSE
          RUN=1    SUBR,(ARGS)
.1        ENDIF
CALL      ENDM
FTN=1     SPACE    2,10
**        FTN=1    SUBR,ARGS       PROCESS FTN FORTRAN ARGUMENTS.
*
*         ENTRY:   SUBR = SUBROUTINE NAME.
*                  ARGS = ARGUMENT LIST.
*
FTN=1     MACRO    SUBR,ARGS
          LOCAL    E,I,J,P,ARGLIST
I         DECMIC   0
          SA0      A1
.1        IFC      NE,$ARGS$$
J         DECMIC   6
          USE      FTN.ARG
ARGLIST   BSS      0
          IRP      ARGS
I         DECMIC   'I'+1
P         ARG=2    (ARGS)
.2        IF       -REG,'P'
          IFC      EQ,$'P'$$,1
P         MICRO    1,,$0LNULL+*+400000B$
          CON      'P'
.2        ELSE
          BSS      1
          USE      *
          R=       X'J','P'
          SA'J'    ARGLIST+'I'-1
          USE      FTN.ARG
J         DECMIC   13D-'J'
.2        ENDIF
          IRP
          CON      0
          USE      *
          SA1      ARGLIST
.1        ENDIF
.1        IF       -MIC,ENTRY
E         MICRO    1,,$0$
.1        ELSE
E         MICRO    1,,$'ENTRY'-2$
.1        ENDIF
+         RJ       =X#SUBR
-         VFD      12/0,18/'E'
          SA1      A0
FTN=1     ENDM
RUN=1     SPACE    2,10
**        RUN=1    SUBR,ARGS       PROCESS RUN FORTRAN ARGUMENTS.
*
*         ENTRY:   SUBR = NAME OF THE SUBROUTINE.
*                  ARGS = LIST OF ARGUMENTS.
*
RUN=1     MACRO    SUBR,ARGS
          LOCAL    E,I,J,P
I         MICRO    1,,$0$
.1        IFC      NE,$ARGS$$
          IRP      ARGS
I         DECMIC   'I'+1
P         ARG=2    (ARGS)
.2        IFLE     'I',6
          IFC      NE,$'P'$B'I'$,2
          IFC      NE,$'P'$$,1
          SB'I'    'P'
.2        ELSE
.3        IFEQ     'I',7
          SA1      =X#SUBR-1
          SB7      X1-6
          SB7      A1-B7
.4        IFC      NE,$'P'$$
.4        IFC      NE,$'P'$X6$
          SX6      'P'
          SA6      B7
.4        ENDIF
J         DECMIC   6
.3        ELSE
.4        IFC      EQ,$'P'$$
          SB7      B7+1
.4        ELSE
J         DECMIC   13D-'J'
          SX'J'    'P'
          SA'J'    B7+1
          SB7      A'J'
.4        ENDIF
.3        ENDIF
.2        ENDIF
          IRP
.1        ENDIF
.1        IF       -MIC,ENTRY
E         MICRO    1,,$0$
.1        ELSE
E         MICRO    1,,$'ENTRY'-1$
.1        ENDIF
+         RJ       =X#SUBR
-         VFD      6/7,6/'I'D,18/'E'
RUN=1     ENDM
ARG=2     SPACE    2,10
** MIC    ARG=2    ARG1            PROCESS CALL MACRO ARGUMENT.
*
*         ENTRY:   MIC = NAME OF THE MICRO TO BE SET TO THE ARGUMENT
*                        STRING ADJUSTED FOR LITERAL VALUES.
*                  ARG1 = ARGUMENT STRING.
*
          MACRO    ARG=2,MIC,ARG1
          LOCAL    C
C         MICRO    1,1,$ARG1$
.1        IFC      NE,$'C'$=$
.1        IFC      GE,$'C'$0$
.2        IFC      LT,$'C'$+$
MIC       MICRO    1,,$=ARG1$
.2        ELSE
.1        IFC      LT,$'C'$*$
C         MICRO    2,,$ARG1$
C         ARG=2    'C'
C         MICRO    1,1,$'C'$
.1        IFC      EQ,$'C'$=$
MIC       MICRO    1,,$=ARG1$
.1        ELSE
MIC       MICRO    1,,$ARG1$
.2        ENDIF
.1        ENDIF
ARG=2     ENDM
          END
*DECK,SDOT
          IDENT  SDOT
*
***       REAL FUNCTION  SDOT(N,SX,INCX,SY,INCY)
*
*         COMPUTED AS SUM FROM I=1 TO N OF  SXII *SYII
*
*         SXII  = SX(1 + (I-1)*INCX)  IF INCX .GE. 0
*               = SX(1 + (I-N)*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR SYII
*
*         SX( ),SY( )               SINGLE PRECISION
*         N,INCX,INCY               INTEGER TYPE
*         SUM ACCUMULATED IN        SINGLE PRECISION
*         RESULT  SDOT  IN          SINGLE PRECISION
*
*         ROUNDED ARITHMETIC INSTRUCTIONS ARE USED
*
*         WRITTEN BY  CLEVE B. MOLER
*                     UNIVERSITY OF NEW MEXICO
*                     ALBUQUERQUE, NEW MEXICO
C
***       1 JUNE 77
*
          ENTRY  SDOT
          VFD    42/4HSDOT,18/5
*
 SDOT     DATA   0
          INFTN  SDOT,5      PROPER LINKAGE (RUN,FTN) MACRO.
*
          MX6    0           (X6)=SDOT=0.
          SA1    B1          (X1)=N
          SB7    1           (B7)=1
*
          SB1    X1          (B1)=N
          SB1    B1-B7       (B1)=N-1
          MI     B1,OUT      IF (N .LE. 0), QUIT
*
          SA1    B2          (X1)=SX(1)
          SA3    B3          (X3)=INCX
*
          SA2    B4          (X2)=SY(1)
          SA4    B5          (X4)=INCY
*
          NZ     B1,NGT1     IF (N .GT. 1), LOOP NEEDED
          RX6    X1*X2       (X6)=SX(1)*SY(1)
          NX6    X6          (X7)=NORM.(X6)
          JP     OUT
 NGT1     SX0    -B1         (X0)=-(N-1)
*
          SB3    X3          (B3)=INCX
          SB4    X4          (B4)=INCY
*
          GE     B3,INCXNN   IF (INCX .GE. 0) NO ADDRESS FIXUP NEEDED
          DX3    X0*X3       (X3)=-(N-1)*INCX
          SB7    A1          (B7)=LOC(SX(1))
          SA1    B7+X3       (X1)=SX(1+(1-N)*INCX). (A1)=LOC(X(1))
*
 INCXNN   SA3    A1+B3       (X3)=SX(2)
          GE     B4,INCYNN   IF (INCY .GE. 0) NO ADDRESS FIXUP NEEDED
          DX4    X0*X4       (X4)=-(N-1)*INCY
          SB7    A2          (B7)=LOC(SY(1))
          SA2    B7+X4       (X2)=SY(1+(1-N)*INCY). (A2)=LOC(Y(1))
 INCYNN   SA4    A2+B4       (X4)=SY(2)
          SB5    1           (B5)=I=1
          SB6    4           (B6)=4
          SB1    B1-B6       (B1)=N-5
*
          MX0    0           (X0)=0.
          MX5    0           (X5)=0.
          MX7    0           (X7)=0.
          GT     B5,B1,CLEAN IF (I .GT. N-5) CLEAN-UP LOGIC
 LOOP     RX6    X1*X2       (X6)=SX(I)*SY(I)
          SA1    A3+B3       (X1)=SX(I+2)
          SA2    A4+B4       (X2)=SY(I+2)
          NX5    X5          (X5)=NORM.(X5)
          RX0    X0+X7       (X0)=SUM1=SUM1+SX(I-1)*SY(I-1)
*
          RX7    X3*X4       (X7)=SX(I+1)*SY(I+1)
          SA3    A1+B3       (X3)=SX(I+3)
          SA4    A2+B4       (X4)=SX(I+3)
          NX0    X0          (X0)=NORM.(X0)
          RX5    X5+X6       (X5)=SUM2=SUM2+SX(I)*SY(I)
*
          SB5    B5+B6       (B5)=I=I+4. INCREMENT I.
          RX6    X1*X2       (X6)=SX(I-2)*SY(I-2)
          SA1    A3+B3       (X1)=SX(I)
          SA2    A4+B4       (X2)=SY(I)
          NX5    X5          (X5)=NORM.(X5)
          RX0    X0+X7       (X0)=SUM1+SX(I-3)*SY(I-3)
*
          RX7    X3*X4       (X7)=SX(I-1)*SY(I-1)
          SA3    A1+B3       (X3)=SX(I+1)
          SA4    A2+B4       (S4)=SY(I+1)
          NX0    X0          (X0)=NORM.(X0)
          RX5    X5+X6       (X5)=SUM2=SUM2+SX(I-2)*SY(I-2)
*
          LE     B5,B1,LOOP  IF (I .LE. N-5) CONTINUE LOOP
 CLEAN    SB6    2           (B6)=2
          SB1    B1+B6       (B1)=N-3
          GT     B5,B1,SWAB  IF (I .GT. N-3) 3 OR LESS COMPS. REMAIN
          RX6    X1*X2       (X6)=SX(I)*SY(I)
          SA1    A3+B3       (X1)=SX(I+2)
          SA2    A4+B4       (X2)=SY(I+2)
          NX5    X5          (X5)=NORM.(X5)
          RX0    X0+X7       (X0)=SUM1=SUM1+SX(I-1)*SY(I-1)
*
          RX7    X3*X4       (X7)=SX(I+1)*SY(I+1)
          SA3    A1+B3       (X3)=SX(I+3)
          SA4    A2+B4       (X4)=SY(I+3)
          RX5    X5+X6       (X5)=SUM2=SUM2+SX(I)*SY(I)
          NX0    X0          (X0)=NORM.(X0)
*
          SB5    B5+B6       (B5)=I=I+2. INCREMENT I
 SWAB     SB1    B1+B6       (B1)=N-1
          GT     B5,B1,MOP   IF (I .GT. N-1) AT MOST 1 COMP. REMAINS
          RX6    X1*X2       (X6)=SX(I)*SY(I)
          NX5    X5          (X5)=NORM.(X5)
          RX0    X0+X7       (X0)=SUM1=SUM1+SX(I-1)*SY(I-1)
*
          RX7    X3*X4       (X7)=SX(I+1)*SY(I+1)
          RX5    X5+X6       (X5)=SUM2=SUM2+SX(I)*SY(I)
          NX0    X0          (X0)=NORM.(X0)
          SB5    B5+B6       (B5)=I=I+2. INCREMENT I.
*
 MOP      SB1    B1+B6       (B1)=N+1
          GE     B5,B1,WIPE  IF (I .GT. N) PUT ODD-EVEN PARTS TOGETHER
          SA1    A3+B3       (X1)=SX(N)
          SA2    A4+B4       (X2)=SY(N)
          RX6    X1*X2       (X6)=SX(N)*SY(N)
          NX5    X5          (X5)=NORM.(X5)
          RX5    X5+X6       (X5)=SUM2=SUM2+SX(N)*SY(N)
 WIPE     RX0    X0+X7       SUM EVEN INDEXED PRODUCTS.
          RX6    X0+X5       (X6)=SUM(SX(I)*SY(I))
          NX6    X6          (X6)=NORM.(X6)
 OUT      OUTFTN SDOT        RETURN
*         END    SDOT
          END
*DECK,DSDOT
          IDENT  DSDOT
*
***       DOUBLE FUNCTION  DSDOT(N,SX,INCX,SY,INCY)
*
*         COMPUTED AS SUM FROM I=1 TO N OF  SXII *SYII
*
*         SXII  = SX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*               = SX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR SYII
*
*         SX( ),SY( )               SINGLE PRECISION
*         N,INCX,INCY               INTEGER TYPE
*         SUM ACCUMULATED IN        DOUBLE PRECISION
*         RESULT  DSDOT  IN         DOUBLE PRECISION
*
*         WRITTEN BY  DAVID R. KINCAID
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  DSDOT
          VFD    42/5HDSDOT,18/5
*
 DSDOT    DATA   0               ENTRY/EXIT
          INFTN  DSDOT,5
          SA1    B1              (X1) = N
          SB7    -1              (B7) = -1
          MX6    0
          SB1    X1+B7           (B1) = N-1
          MX7    0               (X6,X7) = 0
*
          SA3    B3              (X3) = INCX
          NG     B1,OUT          IF N .LE. 0 , GO TO OUT
          SA5    B5              (X5) = INCY
          SX1    -B1             (X1) = -(N-1)
          SB3    X3              (B3) = INCX
          SB5    X5              (B5) = INCY
*
          GT     B3,ONE          IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3           LOC(SXI1 ) = LOC(SX) - (N-1)*INCX
          SB2    X3+B2           (B2) = LOC(SXI1 )
*
 ONE      GT     B5,TWO          IF INCY .GT. 0 , GO TO TWO
          DX5    X1*X5           LOC(SYI1 ) = LOC(SY) - (N-1)*INCY
          SB4    X5+B4           (B4) = LOC(SYI1 )
*
*                                (I=1)
 TWO      SA1    B2              (X1) = SXI1
          SA3    B4              (X3) = SYI1
*
          FX0    X1*X3           (X0,X2) = SXI1 *SYI1
          DX2    X1*X3
*
          ZR     B1,EXIT         IF I .EQ. N , GO TO EXIT
*
*                                (I = I+1)
 LOOP     SA1    A1+B3           (X1) = SXII
          SA3    A3+B5           (X3) = SYII
*
          FX4    X6+X0           (X6,X7) = (X6,X7) + (X0,X2)
          DX5    X6+X0
          FX0    X7+X2
          NX4    X4
          FX2    X0+X5
          FX0    X2+X4
          NX5    X0
          DX2    X2+X4
          NX4    X2
          FX6    X4+X5
          DX7    X4+X5
*
          FX0    X1*X3
          SB1    B1+B7           COUNT TERM
          DX2    X1*X3           (X0,X2) = SXII *SYII
*
          NZ     B1,LOOP         IF I .NE. N , GO TO LOOP
*
*                                (I=N)
 EXIT     FX4    X6+X0           (X6,X7) = (X6,X7) + (X0,X2)
          DX5    X6+X0
          FX0    X7+X2
          NX4    X4
          FX2    X0+X5
          FX0    X2+X4
          NX5    X0
          DX2    X2+X4
          NX4    X2
          FX6    X4+X5
          DX7    X4+X5
*
 OUT      OUTFTN DSDOT           RETURN
          END
*DECK,SDSDOT
          IDENT  SDSDOT
*
***       REAL FUNCTION  SDSDOT(N,SB,SX,INCX,SY,INCY)
*
*         COMPUTED AS SUM FROM I=1 TO N OF  SXII *SYII
*
*         SXII  = SX(1 + (I-1)*INCX)  IF INCX .GE. 0
*               = SX(1 + (I-N)*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR SYII
*
*         SX( ),SY( )               SINGLE PRECISION
*         N,INCX,INCY               INTEGER TYPE
*         SUM ACCUMULATED IN        DOUBLE PRECISION
*         RESULT  SDSDOT IN         SINGLE PRECISION (ROUNDED)
*
*
*         WRITTEN BY  DAVID R. KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  SDSDOT
          VFD    42/6HSDSDOT,18/6
*
 SDSDOT   DATA   0               ENTRY/EXIT
          INFTN  SDSDOT,6
*
          SX6    B2
          SA6    ADRSB         SAVE ADDRESS OF SB
*
          CALL   DSDOT,(B1,B3,B4,B5,B6)
*
          SA4    ADRSB         (X4) = SB
          SA4    X4
          FX1    X4+X6
          DX2    X4+X6
          FX3    X2+X7
          FX2    X1+X3
          NX0    X2
          DX3    X1+X3
          NX1    X3
          FX2    X0+X1
          DX3    X0+X1
          RX2    X2+X3
          NX6    X2
*
 OUT      OUTFTN SDSDOT          RETURN
*
 ADRSB    BSS    1             ADDRESS OF SB
*
          END
*DECK,DDOT
          IDENT  DDOT
*
***       DOUBLE FUNCTION  DDOT(N,DX,INCX,DY,INCY)
*
*         COMPUTED AS SUM FROM I=1 TO N OF  DXII *DYII
*
*         DXII  = DX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*               = DX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR  DYII
*
*         DX( ),DY( )               DOUBLE PRECISION
*         N,INCX,INCY               INTEGER TYPE
*         SUM ACCUMULATED IN        DOUBLE PRECISION
*         RESULT  DDOT  IN          DOUBLE PRECISION
*
*         WRITTEN BY  DAVID R. KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  DDOT
          VFD    42/4HDDOT,18/5
*
 DDOT     DATA   0             ENTRY/EXIT
          INFTN  DDOT,5
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          MX6    0
          SB1    X1+B7         (B1) = N-1
          MX7    0             (X6,X7) = 0
*
          SA3    B3            (X3) = INCX
          NG     B1,OUT        IF N .LE. 0 , GO TO OUT
          SA5    B5            (X5) = INCY
          SX1    -B1           (X1) = -(N-1)
          LX3    1             INCX = 2*INCX
          IX5    X5+X5         INCY = 2*INCY
          SB3    X3            (B3) = INCX
          SB5    X5            (B5) = INCY
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3         LOC(DXI1 ) = LOC(DX) - (N-1)*INCX
          SB2    X3+B2         (B2) = LOC(DXI1 )
*
 ONE      GT     B5,TWO        IF INCY .GT. 0 , GO TO TWO
          DX5    X1*X5         LOC(DYI1 ) = LOC(DY) - (N-1)*INCY
          SB4    X5+B4         (B4) = LOC(DYI1 )
*
*                              (I=1)
 TWO      SA1    B2
          SA3    B4
          SA2    B2-B7         (X1,X2) = DXI1
          SA4    B4-B7         (X3,X4) = DYI1
*
          FX5    X2*X3         (X0,X2) = DXI1 *DYI1
          FX0    X1*X4
          FX5    X0+X5
          FX4    X1*X3
          DX0    X1*X3
          FX5    X0+X5
          FX0    X4+X5
          DX2    X4+X5
*
          ZR     B1,EXIT       IF I .EQ. N , GO TO EXIT
*
*                              (I = I+1)
 LOOP     SA1    A1+B3
          SA3    A3+B5
*
          FX4    X6+X0         (X6,X7) = (X6,X7) + (X0,X2)
          DX5    X6+X0
          FX0    X7+X2
          NX4    X4
          FX2    X0+X5
          FX0    X2+X4
          NX5    X0
          DX2    X2+X4
          NX4    X2
          FX6    X4+X5
          DX7    X4+X5
*
          SB1    B1+B7         COUNT TERM
          SA2    A1-B7         (X1,X2) = DXII
          SA4    A3-B7         (X3,X4) = DYII
*
          FX5    X2*X3         (X0,X2) = DXII *DYII
          FX0    X1*X4
          FX5    X0+X5
          FX4    X1*X3
          DX0    X1*X3
          FX5    X0+X5
          FX0    X4+X5
          DX2    X4+X5
*
          NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
*                              (I=N)
 EXIT     FX4    X6+X0         (X6,X7) = (X6,X7) + (X0,X2)
          DX5    X6+X0
          FX0    X7+X2
          NX4    X4
          FX2    X0+X5
          FX0    X2+X4
          NX5    X0
          DX2    X2+X4
          NX4    X2
          FX6    X4+X5
          DX7    X4+X5
*
 OUT      OUTFTN DDOT          RETURN
          END
*DECK,DQDOTI
          IDENT  DQDOTI
          ENTRY  DQDOTI
 ARG7     BSS    1
          VFD    42/6HDQDOTI,18/7
 DQDOTI   DATA   0
          EQ     DQDOTI
          END
*DECK,DQDOTA
          IDENT  DQDOTA
          ENTRY  DQDOTA
 ARG7     BSS    1
          VFD    42/6HDQDOTA,18/7
 DQDOTA   DATA   0
          EQ     DQDOTA
          END
*DECK,CDOTC
          IDENT  CDOTC
*
***       COMPLEX FUNCTION  CDOTC(N,CX,INCX,CY,INCY)
*
*         COMPUTED AS SUM FROM I=1 TO N OF  CONJ(CXII )*CYII
*
*         CXII  = CX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*               = CX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR  CYII
*
*         CX( ),CY( )               COMPLEX TYPE
*         N,INCX,INCY               INTEGER TYPE
*         SUM ACCUMULATED IN        SINGLE PRECISION
*         RESULT  CDOTC IN          COMPLEX TYPE
*
*         ROUNDED ARITHMETRIC INSTRUCTIONS ARE USED
*
*         WRITTEN BY  DAVID R. KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  CDOTC
          VFD    42/5HCDOTC,18/5
*
 CDOTC    DATA   0             ENTRY/EXIT
          INFTN  CDOTC,5
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          MX6    0
          SB1    X1+B7         (B1) = N-1
          MX7    0             (X6,X7) = 0
*
          SA3    B3            (X3) = INCX
          NG     B1,OUT        IF N .LE. 0 , GO TO OUT
          SA5    B5            (X5) = INCY
          SX1    -B1           (X1) = -(N-1)
          LX3    1             INCX = 2*INCX
          IX5    X5+X5         INCY = 2*INCY
          SB3    X3            (B3) = INCX
          SB5    X5            (B5) = INCY
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3         LOC(CXI1 ) = LOC(CX) - (N-1)*INCX
          SB2    X3+B2         (B2) = LOC(CXI1 )
*
 ONE      GT     B5,TWO        IF INCY .GT. 0 , GO TO TWO
          DX5    X1*X5         LOC(CYI1 ) = LOC(CY) - (N-1)*INCY
          SB4    X5+B4         (B4) = LOC(CYI1 )
*
*                              (I=1)
 TWO      SA1    B2            (X1) = REAL(CXI1 )
          SA2    B4            (X2) = REAL(CYI1 )
*
          RX0    X1*X2         (X0) = REAL(CXI1 )*REAL(CYI1 )
*
          RX5    X6+X0         (X6) = (X6) + (X0)
          NX6    X5
*
          SA4    A2-B7         (X4) = IMAG(CYI1 )
*
          RX0    X1*X4         (X0) = REAL(CXI1 )*IMAG(CYI1 )
*
          RX5    X7+X0         (X7) = (X7) + (X0)
          NX7    X5
*
          SA3    A1-B7         (X3) = IMAG(CXI1 )
*
          RX0    X3*X4         (X0) = IMAG(CXI1 )*IMAG(CYI1 )
*
          RX5    X6+X0         (X6) = (X6) + (X0)
          NX6    X5                 = REAL(CONJ(CXI1 )*CYI1 )
*
          RX0    X3*X2         (X0) = IMAG(CXI1 )*REAL(CYI1 )
*
          RX5    X7-X0         (X7) = (X7) - (X0)
          NX7    X5                 = IMAG(CONJ(CXI1 )*CYI1 )
*
          ZR     B1,OUT        IF I .EQ. N , GO TO OUT
*
*                              (I = I+1)
 LOOP     SA1    A1+B3         (X1) = REAL(CXII )
          SA2    A2+B5         (X2) = REAL(CYII )
*
          RX0    X1*X2         (X0) = REAL(CXII )*REAL(CYII )
*
          RX5    X6+X0         (X6) = (X6) + (X0)
          NX6    X5
*
          SA4    A2-B7         (X4) = IMAG(CYII )
*
          RX0    X1*X4         (X0) = REAL(CXII )*IMAG(CYII )
*
          RX5    X7+X0         (X7) = (X7) + (X0)
          NX7    X5
*
          SA3    A1-B7         (X3) = IMAG(CXII )
*
          RX0    X3*X4         (X0) = IMAG(CXII )*IMAG(CYII )
*
          RX5    X6+X0         (X6) = (X6) + (X0)
          NX6    X5                 = REAL(CONJ(CXII )*CYII )
*
          RX0    X3*X2         (X0) = IMAG(CXII )*REAL(CYII )
          SB1    B1+B7         COUNT TERM
*
          RX5    X7-X0         (X7) = (X7) - (X0)
          NX7    X5                 = IMAG(CONJ(CXII )*CYII )
*
          NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
 OUT      OUTFTN CDOTC         RETURN
          END
*DECK,CDOTU
          IDENT  CDOTU
*
***       COMPLEX FUNCTION  CDOTU(N,CX,INCX,CY,INCY)
*
*         COMPUTED AS SUM FROM I=1 TO N OF  CXII *CYII
*
*         CXII  = CX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*               = CX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR  CYII
*
*         CX( ),CY( )               COMPLEX TYPE
*         N,INCX,INCY               INTEGER TYPE
*         SUM ACCUMULATED IN        SINGLE PRECISION
*         RESULT  CDOTU  IN         COMPLEX TYPE
*
*         ROUNDED ARITHMETRIC INSTRUCTIONS ARE USED
*
*         WRITTEN BY  DAVID R. KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       15 OCT 1974
***       1 JUNE 77
*
          ENTRY  CDOTU
          VFD    42/5HCDOTU,18/5
*
 CDOTU    DATA   0             ENTRY/EXIT
          INFTN  CDOTU,5
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          MX6    0
          SB1    X1+B7         (B1) = N-1
          MX7    0             (X6,X7) = 0
*
          SA3    B3            (X3) = INCX
          NG     B1,OUT        IF N .LE. 0 , GO TO OUT
          SA5    B5            (X5) = INCY
          SX1    -B1           (X1) = -(N-1)
          LX3    1             INCX = 2*INCX
          IX5    X5+X5         INCY = 2*INCY
          SB3    X3            (B3) = INCX
          SB5    X5            (B5) = INCY
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3         LOC(CXI1 ) = LOC(CX) - (N-1)*INCX
          SB2    X3+B2         (B2) = LOC(CXI1 )
*
 ONE      GT     B5,TWO        IF INCY .GT. 0 , GO TO TWO
          DX5    X1*X5         LOC(CYI1 ) = LOC(CY) - (N-1)*INCY
          SB4    X5+B4         (B4) = LOC(CYI1 )
*
*                              (I=1)
 TWO      SA1    B2            (X1) = REAL(CXI1 )
          SA2    B4            (X2) = REAL(CYI1 )
*
          RX0    X1*X2         (X0) = REAL(CXI1 )*REAL(CYI1 )
*
          RX5    X6+X0         (X6) = (X6) + (X0)
          NX6    X5
*
          SA4    A2-B7         (X4) = IMAG(CYI1 )
*
          RX0    X1*X4         (X0) = REAL(CXI1 )*IMAG(CYI1 )
*
          RX5    X7+X0         (X7) = (X7) + (X0)
          NX7    X5
*
          SA3    A1-B7         (X3) = IMAG(CXI1 )
*
          RX0    X3*X4         (X0) = IMAG(CXI1 )*IMAG(CYI1 )
*
          RX5    X6-X0         (X6) = (X6) - (X0)
          NX6    X5                 = REAL(CXI1 *CYI1 )
*
          RX0    X3*X2         (X0) = IMAG(CXI1 )*REAL(CYI1 )
*
          RX5    X7+X0         (X7) = (X7) + (X0)
          NX7    X5                 = IMAG(CXI1 *CYI1 )
*
          ZR     B1,OUT        IF I .EQ. N , GO TO OUT
*
*                              (I = I+1)
 LOOP     SA1    A1+B3         (X1) = REAL(CXII )
          SA2    A2+B5         (X2) = REAL(CYII )
*
          RX0    X1*X2         (X0) = REAL(CXII )*REAL(CYII )
*
          RX5    X6+X0         (X6) = (X6) + (X0)
          NX6    X5
*
          SA4    A2-B7         (X4) = IMAG(CYII )
*
          RX0    X1*X4         (X0) = REAL(CXII )*IMAG(CYII )
*
          RX5    X7+X0         (X7) = (X7) + (X0)
          NX7    X5
*
          SA3    A1-B7         (X3) = IMAG(CXII )
*
          RX0    X3*X4         (X0) = IMAG(CXII )*IMAG(CYII )
*
          RX5    X6-X0         (X6) = (X6) - (X0)
          NX6    X5                 = REAL(CXII *CYII )
*
          RX0    X3*X2         (X0) = IMAG(CXII )*REAL(CYII )
          SB1    B1+B7         COUNT TERM
*
          RX5    X7+X0         (X7) = (X7) + (X0)
          NX7    X5                 = IMAG(CXII *CYII )
*
          NZ     B1,LOOP       IF I .NE. 0 , GO TO LOOP
*
 OUT      OUTFTN CDOTU         RETURN
          END
*DECK,CZDOTC
          IDENT  CZDOTC
*
***       COMPLEX FUNCTION  CZDOTC(N,CX,INCX,CY,INCY)
*
*         COMPUTED AS SUM FROM I=1 TO N OF  CONJ(CXII )*CYII
*
*         CXII  = CX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*               = CX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR  CYII
*
*         CX( ),CY( )               COMPLEX TYPE
*         N,INCX,INCY               INTEGER TYPE
*         SUM ACCUMULATED IN        DOUBLE PRECISION
*         RESULT  CZDOTC IN         COMPLEX TYPE  (ROUNDED)
*
*         WRITTEN BY  DAVID R. KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  CZDOTC
          VFD    42/6HCZDOTC,18/6
*
 CZDOTC   DATA   0             ENTRY/EXIT
          INFTN  CZDOTC,6
          SA1    B1            (X1) = N
          MX4    0
          SB7    -1            (B7) = -1
          SB1    X1            (B1) = N    (I=0)
          MX5    0
          SB6    X1+B7         (B6) = N-1
          BX6    X4
          BX7    X5            (X7,X5) = (X6,X4) = (0,0)
*
          SA3    B3            (X3) = INCX
          NG     B6,OUT        IF N .LE. 0 , GO TO OUT
          SA2    B5            (X2) = INCY
          SX1    -B6           (X1) = -(N-1)
          LX3    1             INCX = 2*INCY
          IX2    X2+X2         INCY = 2*INCY
          SB3    X3            (B3) = INCX
          SB5    X2            (B5) = INCY
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3         LOC(CXI1 ) = LOC(CX) - (N-1)*INCX
          SB2    X3+B2         (B2) = LOC(CXI1 )
*
 ONE      GT     B5,LOOP       IF INCY .GT. 0 , GO TO LOOP
          DX2    X1*X2         LOC(CYI1 ) = LOC(CY) - (N-1)*INCY
          SB4    X2+B4         (B4) = LOC(CXI1 )
*
*                              (I = I+1)
 LOOP     SA1    B2            (X1) = REAL(CXII )
          SB2    B2+B3
          SA2    B4            (X2) = REAL(CYII )
          SB4    B4+B5
*
          FX0    X1*X2         (X0,X1) = REAL(CXII )*REAL(CYII )
          DX1    X1*X2
*
          FX2    X6+X0         (X6,X4) = (X6,X4) + (X0,X1)
          DX3    X6+X0
          FX0    X4+X1
          NX2    X2
          FX1    X0+X3
          FX0    X1+X2
          NX3    X0
          DX1    X1+X2
          NX2    X1
          FX6    X2+X3
          DX4    X2+X3
*
          SA1    A1-B7         (X1) = IMAG(CXII )
          SA2    A2-B7         (X2) = IMAG(CYII )
*
          FX0    X1*X2         (X0,X1) = IMAG(CXII )*IMAG(CYII )
          DX1    X1*X2
*
          FX2    X6+X0         (X6,X4) = (X6,X4) + (X0,X1)
          DX3    X6+X0
          FX0    X4+X1                 = REAL(CONJ(CXII )*CYII )
          NX2    X2
          FX1    X0+X3
          FX0    X1+X2
          NX3    X0
          DX1    X1+X2
          NX2    X1
          FX6    X2+X3
          DX4    X2+X3
*
          SA1    A1            (X1) = IMAG(CXII )
          SA2    A2+B7         (X2) = REAL(CYII )
*
          FX0    X1*X2         (X0,X1) = IMAG(CXII )*REAL(CYII )
          DX1    X1*X2
*
          FX2    X7-X0         (X7,X5) = (X7,X5) - (X0,X1)
          DX3    X7-X0
          FX0    X5-X1
          NX2    X2
          FX1    X0+X3
          FX0    X1+X2
          NX3    X0
          DX1    X1+X2
          NX2    X1
          FX7    X2+X3
          DX5    X2+X3
*
          SA1    A1+B7         (X1) = REAL(CXII )
          SA2    A2-B7         (X2) = IMAG(CYII )
*
          FX0    X1*X2         (X0,X1) = REAL(CXII )*IMAG(CYII )
          DX1    X1*X2
          SB1    B1+B7         COUNT TERM
*
          FX2    X7+X0         (X7,X5) = (X7,X5) + (X0,X1)
          DX3    X7+X0
          FX0    X5+X1                 = IMAG(CONJ(CXII )*CYII )
          NX2    X2
          FX1    X0+X3
          FX0    X1+X2
          NX3    X0
          DX1    X1+X2
          NX2    X1
          FX7    X2+X3
          DX5    X2+X3
*
          NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
          RX0    X6+X4         ROUNDED FINAL RESULT
          RX1    X7+X5
          NX6    X0
          NX7    X1
*
 OUT      OUTFTN CZDOTC        RETURN
          END
*DECK,CZDOTU
          IDENT  CZDOTU
*
***       COMPLEX FUNCTION  CZDOTU(N,CX,INCX,CY,INCY)
*
*         COMPUTED AS SUM FROM I=1 TO N OF  CXII *CYII
*
*         CXII  = CX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*               = CX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR  CYII
*
*         CX( ),CY( )               COMPLEX TYPE
*         N,INCX,INCY               INTEGER TYPE
*         SUM ACCUMULATED IN        DOUBLE PRECISION
*         RESULT  CZDOTU  IN        COMPLEX TYPE  (ROUNDED)
*
*         WRITTEN BY  DAVID R. KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  CZDOTU
          VFD    42/6HCZDOTU,18/5
*
 CZDOTU   DATA   0             ENTRY/EXIT
          INFTN  CZDOTU,5
          SA1    B1            (X1) = N
          MX4    0
          SB7    -1            (B7) = -1
          SB1    X1            (B1) = N    (I=0)
          MX5    0
          SB6    X1+B7         (B6) = N-1
          BX6    X4
          BX7    X5            (X7,X5) = (X6,X4) = (0,0)
*
          SA3    B3            (X3) = INCX
          NG     B6,OUT        IF N .LE. 0 , GO TO OUT
          SA2    B5            (X2) = INCY
          SX1    -B6           (X1) = -(N-1)
          LX3    1             INCX = 2*INCY
          IX2    X2+X2         INCY = 2*INCY
          SB3    X3            (B3) = INCX
          SB5    X2            (B5) = INCY
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3         LOC(CXI1 ) = LOC(CX) - (N-1)*INCX
          SB2    X3+B2         (B2) = LOC(CXI1 )
*
 ONE      GT     B5,LOOP       IF INCY .GT. 0 , GO TO LOOP
          DX2    X1*X2         LOC(CYI1 ) = LOC(CY) - (N-1)*INCY
          SB4    X2+B4         (B4) = LOC(CXI1 )
*
*                              (I = I+1)
 LOOP     SA1    B2            (X1) = REAL(CXII )
          SB2    B2+B3
          SA2    B4            (X2) = REAL(CYII )
          SB4    B4+B5
*
          FX0    X1*X2         (X0,X1) = REAL(CXII )*REAL(CYII )
          DX1    X1*X2
*
          FX2    X6+X0         (X6,X4) = (X6,X4) + (X0,X1)
          DX3    X6+X0
          FX0    X4+X1
          NX2    X2
          FX1    X0+X3
          FX0    X1+X2
          NX3    X0
          DX1    X1+X2
          NX2    X1
          FX6    X2+X3
          DX4    X2+X3
*
          SA1    A1-B7         (X1) = IMAG(CXII )
          SA2    A2-B7         (X2) = IMAG(CYII )
*
          FX0    X1*X2         (X0,X1) = IMAG(CXII )*IMAG(CYII )
          DX1    X1*X2
*
          FX2    X6-X0         (X6,X4) = (X6,X4) - (X0,X1)
          DX3    X6-X0
          FX0    X4-X1                 = REAL(CXII *CYII )
          NX2    X2
          FX1    X0+X3
          FX0    X1+X2
          NX3    X0
          DX1    X1+X2
          NX2    X1
          FX6    X2+X3
          DX4    X2+X3
*
          SA1    A1            (X1) = IMAG(CXII )
          SA2    A2+B7         (X2) = REAL(CYII )
*
          FX0    X1*X2         (X0,X1) = IMAG(CXII )*REAL(CYII )
          DX1    X1*X2
*
          FX2    X7+X0         (X7,X5) = (X7,X5) + (X0,X1)
          DX3    X7+X0
          FX0    X5+X1
          NX2    X2
          FX1    X0+X3
          FX0    X1+X2
          NX3    X0
          DX1    X1+X2
          NX2    X1
          FX7    X2+X3
          DX5    X2+X3
*
          SA1    A1+B7         (X1) = REAL(CXII )
          SA2    A2-B7         (X2) = IMAG(CYII )
*
          FX0    X1*X2         (X0,X1) = REAL(CXII )*IMAG(CYII )
          DX1    X1*X2
          SB1    B1+B7         COUNT TERM
*
          FX2    X7+X0         (X7,X5) = (X7,X5) + (X0,X1)
          DX3    X7+X0
          FX0    X5+X1                 = IMAG(CXII *CYII )
          NX2    X2
          FX1    X0+X3
          FX0    X1+X2
          NX3    X0
          DX1    X1+X2
          NX2    X1
          FX7    X2+X3
          DX5    X2+X3
*
          NZ     B1,LOOP       IF I .NE. 0 , GO TO LOOP
*
          RX0    X6+X4         ROUNDED FINAL RESULT
          RX1    X7+X5
          NX6    X0
          NX7    X1
*
 OUT      OUTFTN CZDOTU         RETURN
          END
*DECK,SAXPY
          IDENT  SAXPY
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL SAXPY(N,SA,SX,INCX,SY,INCY)
*
*         SA*SXII  + SYII   REPLACES  SYII   FOR I=1,N
*
*         SXII  = SX(1 + (I-1)*INCX)  IF INCX .GE. 0
*               = SX(1 + (I-N)*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR SYII
*
*         SX( ),SY( )               SINGLE PRECISION
*         N,INCX,INCY               INTEGER TYPE
*         SA                        SINGLE PRECISION
*
*         ROUNDED ARITHMETIC INSTRUCTIONS ARE USED
*
*         WRITTEN BY  RICHARD J. HANSON
*                     SANDIA LABORATORIES
*                     ALBUQUERQUE, NEW MEXICO
***       1 JUNE 77
*
          ENTRY  SAXPY
          VFD    42/5HSAXPY,18/6
*
 SAXPY    DATA   0
          INFTN  SAXPY,6     PROPER LINKAGE (RUN,FTN) MACRO.
          SA1    B1          (X1)=N
          SB7    1           (B7)=1
*
          SB1    X1          (B1)=N
          SB1    B1-B7       (B1)=N-1
          MI     B1,OUT      IF(N .LE. 0), QUIT.
*
          SA5    B2          (X5)=SA
          ZR     X5,OUT      IF(SA .EQ. 0.), QUIT
*
          SA1    B3          (X1)=SX(1)
          SA2    B5          (X2)=SY(1)
          SA3    B4          (X3)=INCX
*
          SA4    B6          (X4)=INCY
*
          NZ     B1,NGT1     IF (N .GT. 1), LOOP NEEDED
          RX6    X1*X5       (X6)=SA*SX(1)
          RX6    X2+X6       (X6)=SA*SX(1)+SY(1)
          NX6    X6          (X6)=NORM.(X6)
          SA6    A2          SY(1)=(X6)
          JP     OUT         QUIT
 NGT1     SX0    -B1         (X0)=-(N-1)
*
          SB3    X3          (B3)=INCX
          SB4    X4          (B4)=INCY
*
          GE     B3,INCXNN   IF (INCX .GE. 0) NO ADDRESS FIXUP NEEDED
          DX3    X0*X3       COMPUTE -(N-1)*INCX
          SB7    A1          (B7)=LOC(SX(1))
          SA1    B7+X3       (X1)=SX(1+(1-N)*INCX). (A1)=LOC(X(1))
*
 INCXNN   SA3    A1+B3       (X3)=SX(2)
          GE     B4,INCYNN   IF (INCY .GE. 0) NO ADDRESS FIXUP NEEDED
          DX4    X0*X4       COMPUTE -(N-1)*INCY
          SB7    A2          (B7)=LOC(SY(1))
          SA2    B7+X4       (X2)=SY(1+(1-N)*INCY). (A2)=LOC(Y(1))
*
 INCYNN   SA4    A2+B4       (X4)=SY(2)
          SB5    1           (B5)=I=1
          SB6    4           (B6)=4
          SA0    A2-B4       (A0)=LOC(Y(1))-INCY
          SB1    B1-B6       (B1)=N-5
*
          GT     B5,B1,CLEAN IF (I .GT. N-5) CLEAN-UP LOGIC
 LOOP     RX6    X1*X5       (X6)=SA*SX(I)
          SA1    A3+B3       (X1)=SX(I+2)
          RX7    X3*X5       (X7)=SA*SX(I+1)
          NO     0           DEAD
          SA3    A1+B3       (X3)=SX(I+3)
          NO     0           DEAD
          RX6    X2+X6       (X6)=SA*SX(I)+SY(I)
          RX7    X4+X7       (X7)=SA*SX(I+1)+SY(I+1)
          SA2    A4+B4       (X2)=SY(I+2)
          RX0    X1*X5       (X0)=SA*SX(I+2)
          NX6    X6          (X6)=NORM.(X6)
          SA4    A2+B4       (X4)=SY(I+3)
          NX7    X7          (X7)=NORM.(X7)
          RX3    X3*X5       (X3)=SA*SX(I+3)
*
          SA1    A3+B3       (X1)=SX(I+4). NEXT ITER.
          SA6    A0+B4       SY(I)=(X6)
          RX6    X0+X2       (X6)=SA*SX(I+2)+SY(I+2)
          SA2    A4+B4       (X2)=SY(I+4). NEXT ITER.
          SA7    A6+B4       SY(I+1)=(X7)
          RX7    X3+X4       (X7)=SA*SX(I+3)+SY(I+3)
          SA3    A1+B3       (X3)=SX(I+5). NEXT ITER.
          NX6    X6          (X6)=NORM.(X6)
          SA4    A2+B4       (X4)=SY(I+5). NEXT ITER.
          NX7    X7          (X7)=NORM.(X7)
          SA6    A7+B4       SY(I+2)=(X6)
          SB5    B5+B6       I=I+4. INCREMENT I
          SA7    A6+B4       SY(I-1)=(X7)
          SA0    A7          ADVANCE ADDRESS OF SY(I+4) FOR NEXT ITER.
          LE     B5,B1,LOOP  IF(I.LE.N-5) CONTINUE LOOP
*
 CLEAN    SB6    2           (B6)=2
          SB1    B1+B6       (B1)=N-3
          GT     B5,B1,SWAB  IF (I .GT. N-3) 3 OR LESS COMPS. REMAIN
          RX6    X1*X5       (X6)=SA*SX(I)
          SA1    A3+B3       (X1)=SX(I+2)
          RX7    X3*X5       (X7)=SA*SX(I+1)
          SA3    A1+B3       (X3)=SX(I+3)
          RX6    X2+X6       (X6)=SA*SX(I)+SY(I)
          RX7    X4+X7       (X7)=SA*SX(I+1)+SY(I+1)
          SA2    A4+B4       (X2)=SY(I+2)
          NX6    X6          (X6)=NORM.(X6)
          SA4    A2+B4       (X4)=SY(I+3)
          NX7    X7          (X7)=NORM.(X7)
*
          SB5    B5+B6       I=I+2. INCREMENT I.
          SA6    A0+B4       SY(I-2)=(X6)
          SA7    A6+B4       SY(I-1)=(X7)
          SA0    A7          ADVANCE ADDRESS TO SY(I)
*
 SWAB     SB1    B1+B6       (B1)=N-1
          GT     B5,B1,MOP   IF (I .GT. N-1) AT MOST 1 COMP. REMAINS
          RX6    X1*X5       (X6)=SA*SX(I)
          RX7    X3*X5       (X7)=SA*SX(I+1)
          SB5    B5+B6       I=I+2. INCREMENT I
          RX6    X2+X6       (X6)=SA*SX(I-2)+SY(I-2)
          RX7    X4+X7       (X7)=SA*SX(I-1)+SY(I-1)
          NX6    X6          (X6)=NORM.(X6)
          NX7    X7          (X7)=NORM.(X7)
          SA6    A0+B4       SY(I-2)=(X6)
          SA7    A6+B4       SY(I-1)=(X7)
          SA0    A7          ADVANCE ADDRESS TO SY(I)
*
 MOP      SB1    B1+B6       (B1)=N+1
          GE     B5,B1,OUT   IF (I .GT. N) RETURN
          SA1    A3+B3       (X1)=SX(N)
          SA2    A4+B4       (X2)=SY(N)
          RX6    X1*X5       (X6)=SA*SX(N)
          RX6    X2+X6       (X6)=SA*SX(N)+SY(N)
          NX6    X6          (X6)=NORM.(X6)
          SA6    A0+B4       SY(N)=(X6)
*
 OUT      OUTFTN SAXPY       RETURN
*         END    SAXPY
          END
*DECK,DAXPY
          IDENT  DAXPY
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL DAXPY(N,DA,DX,INCX,DY,INCY)
*
*         DA*DXII  + DYII   REPLACES  DYII   FOR I=1,N
*
*         DXII  = DX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*               = DX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR DYII
*
*         DX( ),DY( )               DOUBLE PRECISION
*         N,INCX,INCY               INTEGER TYPE
*         DA                        DOUBLE PRECISION
*
*         WRITTEN BY  DAVID R. KINCAID AND ELIZABETH WILLIAMS
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  DAXPY
          VFD    42/5HDAXPY,18/6
*
 DAXPY    DATA   0             ENTRY/EXIT
          INFTN  DAXPY,6
          SA3    B1            (X3) = N
          SB7    -1            (B7) = -1
          SB1    X3+B7         (B1) = N-1
          SA1    B2            (X1,X2) = DA
          SA2    B2-B7
          ZR     X1,OUT        IF(DA .EQ. 0) GO TO OUT
*
          SA4    B4            (X4) = INCX
          NG     B1,OUT        IF N .LE. 0 , GO TO OUT
          SA5    B6            (X5) = INCY
          SX3    -B1           (X3) = -(N-1)
          LX4    1             INCX = 2*INCX
          IX5    X5+X5         INCY = 2*INCY
          SB4    X4            (B4) = INCX
          SB6    X5            (B5) = INCY
*
          GT     B4,ONE        IF INCX .GT. 0 , GO TO ONE
          DX4    X3*X4         LOC(DXI1 ) = LOC(DX) - (N-1)*INCX
          SB3    X4+B3         (B3) = LOC(DXI1 )
*
 ONE      GT     B6,TWO        IF INCY .GT. 0 , GO TO TWO
          DX5    X3*X5         LOC(DYI1 ) = LOC(DY) - (N-1)*INCY
          SB5    X5+B5         (B5) = LOC(DYI1 )
*
*                              (I = 1)
 TWO      SA3    B3            (X3,X4) = DXI1
          SA4    B3-B7
*
          FX5    X2*X3         (X6,X7) = DA*DXI1
          FX0    X1*X4
          FX5    X0+X5
          FX4    X1*X3
          DX0    X1*X3
          FX5    X0+X5
          FX6    X4+X5
          DX7    X4+X5
*
          SA5    B5            (X5,X4) = DYI1
          SA4    B5-B7
*
          FX0    X6+X5         (X6,X7) = (X6,X7) + DYI1
          DX6    X6+X5
          FX5    X7+X4
          NX0    X0
          FX4    X5+X6
          FX5    X4+X0
          NX6    X5
          DX4    X4+X0
          NX0    X4
          FX6    X0+X6
          DX7    X0+X6
*
          SA6    A5            DYI1  = (X6,X7)
          SA7    A4
*
          ZR     B1,OUT        IF I .EQ. N , GO TO OUT
*
*                              (I = I+1)
 LOOP     SA3    A3+B4         (X3,X4) = DXII
          SA4    A3-B7
*
          FX5    X2*X3         (X6,X7) = DA*DXII
          FX0    X1*X4
          FX5    X0+X5
          FX4    X1*X3
          DX0    X1*X3
          FX5    X0+X5
          FX6    X4+X5
          DX7    X4+X5
*
          SA5    A5+B6         (X5,X4) = DYII
          SA4    A5-B7
*
          FX0    X6+X5         (X6,X7) = (X6,X7) + DYII
          DX6    X6+X5
          FX5    X7+X4
          NX0    X0
          FX4    X5+X6
          FX5    X4+X0
          NX6    X5
          DX4    X4+X0
          NX0    X4
          FX6    X0+X6
          DX7    X0+X6
*
          SB1    B1+B7         I = I+1
*
          SA6    A5            DYII  = (X6,X7)
          SA7    A4
*
          NZ     B1,LOOP       IF I .EQ. N , GO TO LOOP
*
 OUT      OUTFTN DAXPY         RETURN
          END
*DECK,CAXPY
          IDENT  CAXPY
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL CAXPY(N,CA,CX,INCX,CY,INCY)
*
*         CA*CXII  + CYII   REPLACES  CYII   FOR I=1,N
*
*         CXII  = CX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*                = CX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR CYII
*
*         CX( ),CY( )               COMPLEX TYPE
*         N,INCX,INCY               INTEGER TYPE
*         CA                        COMPLEX TYPE
*
*         WRITTEN BY  DAVID R. KINCAID AND ELIZABETH WILLIAMS
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  CAXPY
          VFD    42/5HCAXPY,18/6
*
 CAXPY    DATA   0             ENTRY/EXIT
          INFTN  CAXPY,6
          SA3    B1            (X3) = N
          SB7    -1            (B7) = -1
          SB1    X3+B7         (B1) = N-1
*
          SA4    B4            (X4) = INCX
          NG     B1,OUT        IF N .LE. 0 , GO TO OUT
          SA5    B6            (X5) = INCY
          SX3    -B1           (X3) = -(N-1)
          LX4    1             INCX = 2*INCX
          IX5    X5+X5         INCY = 2*INCY
          SB4    X4            (B4) = INCX
          SB6    X5            (B6) = INCY
*
          GT     B4,ONE        IF INCX .GT. 0 , GO TO ONE
          DX4    X3*X4         LOC(CXI1 ) = LOC(CX) - (N-1)*INCX
          SB3    X4+B3         (B3) = LOC(CXI1 )
*
 ONE      GT     B6,TWO        IF INCY .GT. 0 , GO TO TWO
          DX5    X3*X5         LOC(CYI1 ) = LOC(CY) - (N-1)*INCY
          SB5    X5+B5
*
*                              (I = 1)
 TWO      SA3    B3            (X3) = REAL(CXI1 )
          SA1    B2            (X1) = REAL(CA)
          SA2    B2-B7         (X2) = IMAG(CA)
*
          BX5    X1
          AX5    59
          BX5    X1-X5         (X5) = ABS(REAL(CA))
          BX6    X2
          AX6    59
          BX6    X2-X6         (X6) = ABS(IMAG(CA))
          RX6    X5+X6
          NX6    X6
          ZR     X6,OUT        IF(ABS(REAL(CA))+ABS(IMAG(CCA))=0.0) GOTO
*
          SA4    B3-B7         (X4) = IMAG(CXI1 )
*                              (X6,X7) = CA*CXI1
          FX0    X1*X3         (X0) = REAL(CA)*REAL(CXI1 )
          FX5    X2*X4         (X5) = IMAG(CA)*IMAG(CXI1 )
          FX6    X0-X5         (X6) = REAL(CA*CXI1 )
*
          FX0    X1*X4         (X0) = REAL(CA)*IMAG(CXI1 )
          FX5    X2*X3         (X5) = IMAG(CA)*REAL(CXI1 )
          FX7    X0+X5         (X7) = IMAG(CA*CXI1 )
*
*
*                              (X6,X7) = (X6,X7) + CYI1
          SA5    B5            (X5) = REAL(CYI1 )
          SA4    B5-B7         (X4) = IMAG(CYI1 )
*
          FX0    X6+X5         (X0) = REAL(CA*CXI1 ) + REAL(CYI1 )
          FX3    X7+X4         (X3) = IMAG(CA*CXI1 ) + IMAG(CYI1 )
          NX6    X0
          NX7    X3            NORMALIZE RESULT
*
          SA6    A5            REAL(CYI1 ) = (X6)
          SA7    A4            IMAG(CYI1 ) = (X7)
*
          ZR     B1,OUT        IF I .EQ. N , GO TO OUT
*
*                              (I = I+1)
 LOOP     SA3    A3+B4         (X3) = REAL(CXII )
          SA4    A3-B7         (X4) = IMAG(CXII )
*
*                              (X6,X7) = CA*CXII
          FX0    X1*X3         (X0) = REAL(CA)*REAL(CXII )
          FX5    X2*X4         (X5) = IMAG(CA)*IMAG(CXII )
          FX6    X0-X5         (X6) = REAL(CA*CXII )
*
          FX0    X1*X4         (X0) = REAL(CA)*IMAG(CXII )
          FX5    X2*X3         (X5) = IMAG(CA)*REAL(CXII )
          FX7    X0+X5         (X7) = IMAG(CA*CXII )
*
*                              (X6,X7) = (X6,X7) + CYII
          SA5    A5+B6         (X5) = REAL(CYII )
          SA4    A5-B7         (X4) = IMAG(CYII )
*
          FX0    X6+X5         (X0) = REAL(CA*CXII ) + REAL(CYII )
          FX3    X7+X4         (X3) = IMAG(CA*CXII ) + IMAG(CYII )
          SB1    B1+B7         I = I+1
          NX6    X0
          NX7    X3            NORMALIZE RESULT
*
          SA6    A5            REAL(CYII ) = (X6)
          SA7    A4            IMAG(CYII ) = (X7)
*
          NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
 OUT      OUTFTN CAXPY         RETURN
          END
*DECK,SROTG
          IDENT  SROTG
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL SROTG(SA,SB,SC,SS)
*
*         COMPUTE QUANTITIES :     R = SQRT( SA**2 + SB**2 )
*                                  SC = SA/R  ,  SS = SB/R
*                                  SA = R
*
*         DEFINING THE GIVENS REFLECTION MATRIX   (SC   SS)
*                                                 (-SS  SC)
*
*         SA,SB,SC,SS               SINGLE PRECISION
*
*         ROUNDED ARITHMETRIC INSTRUCTIONS ARE USED
*
*
*         WRITTEN BY  DAVID R. KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  SROTG
          VFD    42/5HSROTG,18/4
*
 SROTG    DATA   0             ENTRY/EXIT
          INFTN  SROTG,4
*
          SX6    B1
          SX7    B2
          SA6    ADRSA         SAVE ADDRESS OF SA AND SB
          SA7    ADRSB
          SX6    B3
          SX7    B4
          SA6    ADRSC         SAVE ADDRESS OF SC AND SS
          SA7    ADRSS
          SA2    B1            (X2) = SA
          SA3    B2            (X3) = SB
          SA5    UNIT          (X5) = 1.0
*
          BX4    X3
          AX4    59
          BX7    X3-X4         (X7) = ABS(SB)
          ZR     X7,THIRTY     IF ABS(SB) .EQ. 0 , GO TO THIRTY
          BX1    X2
          AX1    59
          BX6    X2-X1         (X6) = ABS(SA)
          ZR     X6,FORTY      IF ABS(SA) .EQ. 0 , GO TO FORTY
*
          RX6    X6-X7
          NX6    X6
          ZR     X6,TWENTY
          NG     X6,TWENTY     IF ABS(SA) .LE. ABS(SB), GO TO TWENTY
*
          RX6    X3/X2         (X6) = SB/SA   (=XR)
          RX0    X6*X6         (X0) = XR**2
          SA6    XR            XR = (X6)
          RX0    X5+X0
          NX6    X0            (X6) = 1.+XR**2
          SA6    XR2P1         XR2P1 = (X6)
*
          SB1    XR2P1
*
          CALL   SQRT,(B1)     (X6) =   SQRT(1.+XR**2)  (=YR)
*
          SA1    ADRSA         RESTORE B REGISTERS
          SB1    X1
          SA2    ADRSB
          SB2    X2
          SA3    ADRSC
          SB3    X3
          SA4    ADRSS
          SB4    X4
*
          SA2    B1            (X2)=SA
          SA3    XR            (A3) = XR
*
          RX7    X2*X6         (X7) = SA*YR
          SA5    UNIT
          RX6    X5/X6         (X6) = 1./YR
          SA7    B1            SA=(X7)
          RX7    X6*X3         (X7) = SC*XR
          SA6    B3            SC=(X6)
          SA7    B4            SS=(X7)
*
          EQ     FIFTY         GO TO FIFTY
*
 TWENTY   RX7    X2/X3         (X7) = SA/SB  (= XR)
          RX0    X7*X7         (X0) = XR**2
          SA7    XR            XR = (X7)
          RX7    X5+X0
          NX6    X7            (X6) = 1.+XR**2
          SA6    XR2P1         XR2P1 = (X6)
*
          SB1    XR2P1
*
          CALL   SQRT,(B1)     (X6) =   SQRT(1.+XR**2)  (=YR)
*
          SA1    ADRSA         RESTORE B REGISTERS
          SB1    X1
          SA2    ADRSB
          SB2    X2
          SA3    ADRSC
          SB3    X3
          SA4    ADRSS
          SB4    X4
*
          SA3    B2            (X3)=SB
          SA1    XR            (X1) = XR
*
          RX7    X3*X6         (X7) = SB*YR
          SA5    UNIT
          RX6    X5/X6         (X6) = 1./YR
          SA7    B1            SA=(X7)
          RX7    X6*X1         (X7) = SS*XR
          SA6    B4            SS=(X6)
          SA7    B3            SC=(X7)
*
          EQ     FIFTY         GO TO FIFTY
*
 THIRTY   BX6    X5            (X6) = 1.
          MX7    0             (X7) = 0.
          SA6    B3            SC = (X6)
          SA7    B4            SS = (X7)
*
          EQ     FIFTY         GO TO FIFTY
*
 FORTY    BX6    X5            (X6) = 1.
          MX7    0             (X7) = 0.
          SA6    B4            SS = (X6)
          SA7    B3            SC = (X7)
*
          BX6    X3            (X6) = SB
          SA6    B1            SA = (X1)
*
 FIFTY    SA2    B4            (X2) = SS
          SA3    B3            (X3) = SC
          SA5    UNIT          (X5) = 1.0
          ZR     X3,SEVENTY    IF SC .EQ. 0 , TO GO SEVENTY
*
          BX1    X2
          AX1    59
          BX6    X2-X1         (X6) = ABS(SS)
          BX4    X3
          AX4    59
          BX7    X3-X4         (X7) = ABS(SC)
*
          RX6    X6-X7
          NX6    X6
          NG     X6,SIXTY      IF ABS(SS) .LT. ABS(SC), GO TO SIXTY
*
          RX6    X5/X3         (X6) = 1./SC
          SA6    B2            SB = (X6)
          EQ     OUT           GO TO OUT
*
 SIXTY    BX6    X2            (X6) = SC
          SA6    B2            SB = (X6)
          EQ     OUT           GO TO OUT
*
 SEVENTY  BX6    X5            (X6) = 1.
          SA6    B2            SB = (X6)
          EQ     OUT           GO TO OUT
*
 OUT      OUTFTN SROTG         RETURN
*
 ADRSA    BSS    1
 ADRSB    BSS    1
 ADRSC    BSS    1
 ADRSS    BSS    1
*
 XR       BSS    1
 XR2P1    BSS    1
*
 UNIT     DATA   1.0
*
          END
*DECK,DROTG
          IDENT  DROTG
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL DROTG(DA,DB,DC,DS)
*
*         COMPUTE QUANTITIES:   DR = DSQRT( DA**2 + DB**2 )
*                               DC = DA/DR  ,  DS = DB/DR
*                               DA = DR
*
*         DEFINES THE GIVENS REFLECTION MATRIX   (DC   DS)
*                                                (-DS  DC)
*
*         DA,DB,DC,DS                DOUBLE PRECISION
*
*
*         WRITTEN BY  DAVID R. KINCAID AND JAMES SULLIVAN
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  DROTG
          VFD    42/5HDROTG,18/4
 DROTG    DATA   0             ENTRY/EXIT
          INFTN  DROTG,4
*
          SX6    B1
          SX7    B2
          SA6    ADRDA         SAVE ADDRESS OF DA AND DB
          SA7    ADRDB
          SX6    B3
          SX7    B4
          SA6    ADRDC         SAVE ADDRESS OF DC AND DS
          SA7    ADRDS
*
          SB7    -1            (B7) = -1
*
          SA3    B2            (X3,) = DB
          BX7    X3            (X7) = X3
          AX7    73B           FILL X7 WITH THE SIGN BIT OF DB.
          BX4    X7-X3         (X4,) = DABS(DB)
          ZR     X4,THIRTY     IF(SNGL(ABS(DB)) = 0) GO TO THIRTY
*
          SA1    B1            (X1,) = DA
          BX2    X1            (X2) = X1
          BX6    X1            (X6) = X1
          AX6    73B           FILL X6 WITH THE SIGN BIT OF DA.
          BX2    X6-X1         (X2,) = DABS(DA)
*
          ZR     X2,FORTY      IF(SNGL(ABS(DA)) = 0) GO TO FORTY
          FX5    X4-X2         COMPARE UPPER HALVES OF DABS(DA) AND DABS
          NX5    X5            MAKE SURE X5 DOES NOT CONTAIN A MINUS ZER
          NG     X5,TEN        IF (DABS(DA) > DABS(DB)) GO TO TEN.
*                              ELSE IF (SNGL(DABS(DA)) @ SNGL(DABS(DB)))
*                              FOLLOWING....
          SA2    B1-B7         (X1,X2) = DA
          SA4    B2-B7         (X3,X4) = DB
*
          FX5    X1/X3         (X6,X7) = DA / DB
          FX6    X3*X5
          FX7    X1-X6
          DX6    X1-X6
          NX7    X7
          FX6    X6+X7
          DX7    X3*X5
          FX0    X4*X5
          FX6    X2+X6
          FX6    X6-X7
          FX6    X6-X0
          FX0    X6/X3
          FX6    X0+X5
          DX7    X0+X5
          NX5    X6
          FX6    X5+X7
          DX7    X5+X7         (X6,X7) = (X1,X2) / (X3,X4)
*
          SA6    XR            (XR) = (X6,X7)
          SA7    XR+1          (XR) = DA / DB
*
          FX4    X6*X7         (X0,X1) = XR**2
          DX5    X6*X6
          FX4    X4+X4
          FX1    X6*X6
          FX5    X4+X5
          FX0    X1+X5
          DX1    X1+X5         (X0,X1) = (X6,X7) * (X6,X7)
*
          SA4    ONE           (X4) = +1.
*
          FX2    X0+X4         (X6,X7) = 1.D0 + (XR*XR)
          DX3    X0+X4
          NX2    X2
          FX5    X1+X3
          FX4    X2+X5
          NX3    X4
          DX5    X2+X5
          NX2    X5
          FX6    X2+X3
          DX7    X2+X3         (X6,X7) = (1.,0) + (X0,X1)
*
          SA6    XR2P1
          SA7    XR2P1+1
*                              SET (XR2P1) = (X6,X7).
          SB1    XR2P1
          CALL   DSQRT,(B1)
*
          SA1    ADRDA
          SB1    X1            RESTORE B  REGISTERS
          SA2    ADRDB
          SB2    X2
          SA3    ADRDC
          SB3    X3
          SA4    ADRDS
          SB4    X4
          SB7    -1
*                              NO ERROR CHECKS ARE MADE UPON RETURN
          SA6    YR
          SA7    YR+1          (YR) = SGN(B)*DSQRT(ONE+XR*XR)
*
          SA2    ONE           (X2) = +1.
*
          FX1    X2/X6         (X6,X7) = 1.0D0 / YR
          FX4    X1*X6
          FX5    X2-X4
          DX4    X2-X4
          NX5    X5
          FX4    X4+X5
          DX5    X1*X6
          FX0    X1*X7
          FX4    X4-X5
          FX4    X4-X0
          FX0    X4/X6
          FX4    X0+X1
          DX5    X0+X1
          NX1    X4
          FX6    X1+X5
          DX7    X1+X5         (X6,X7) = (X2,) / (X6,X7)
          SA6    B4            DS=(X6,X7)
          SA7    B4-B7
*
          SA2    XR            (X2,X3) = XR
          SA3    XR+1
*
          FX4    X2*X7         (X6,X7) = DS * XR
          FX5    X3*X6
          FX4    X4+X5
          FX7    X2*X6
          DX5    X2*X6
          FX5    X4+X5
          FX6    X5+X7
          DX7    X5+X7         (X6,X7) = (X6,X7) * (X2,X3)
*
          SA6    B3            DC=(X6,X7)
          SA7    B3-B7
*
          SA2    B2            (X2,X3)=DB
          SA3    B2-B7
*
          SA4    YR            (X4,X5) = YR
          SA5    YR+1
*
          FX0    X3*X4         (X6,X7) = DB * YR
          FX1    X2*X5
          FX0    X0+X1
          FX3    X2*X4
          DX1    X2*X4
          FX1    X0+X1
          FX6    X1+X3
          DX7    X1+X3         (X6,X7) = (X2,X3) * (X4,X5)
*
          SA6    B1            DA=(X6,X7)
          SA7    B1-B7
*
          EQ     FIFTY         GO TO FIFTY
*
 TEN      SA2    B1-B7         (X1,X2) = DA
          SA4    B2-B7         (X3,X4) = DB
*
          FX5    X3/X1         (X6,X7) = DB / DA
          FX6    X1*X5
          FX7    X3-X6
          DX6    X3-X6
          NX7    X7
          FX6    X6+X7
          DX7    X1*X5
          FX0    X2*X5
          FX6    X4+X6
          FX6    X6-X7
          FX6    X6-X0
          FX0    X6/X1
          FX6    X0+X5
          DX7    X0+X5
          NX5    X6
          FX6    X5+X7
          DX7    X5+X7         (X6,X7) = (X3,X4) / (X1,X2)
*
          SA6    XR            (XR) = (X6,X7)
          SA7    XR+1          (XR) = DB / DA
*
          FX4    X6*X7         (X0,X1) = XR**2
          DX5    X6*X6
          FX4    X4+X4
          FX3    X6*X6
          FX5    X4+X5
          FX0    X3+X5
          DX1    X3+X5         (X0,X1) = (X6,X7) * (X6,X7)
*
          SA4    ONE           (X4) = +1.
*
          FX2    X0+X4         (X6,X7) = 1.D0 + (XR*XR)
          DX3    X0+X4
          NX2    X2
          FX5    X1+X3
          FX4    X2+X5
          NX3    X4
          DX5    X2+X5
          NX2    X5
          FX6    X2+X3
          DX7    X2+X3         (X6,X7) = (1.,0) + (X0,X1)
*
*
          SA6    XR2P1
          SA7    XR2P1+1
*                              SET (XR2P1) = (X6,X7).
          SB1    XR2P1
          CALL   DSQRT,(B1)
*
          SA1    ADRDA         RESTORE B REGISTERS
          SB1    X1
          SA2    ADRDB
          SB2    X2
          SA3    ADRDC
          SB3    X3
          SA4    ADRDS
          SB4    X4
          SB7    -1
*
*                              NO ERROR CHECKS ARE MADE UPON RETURN
          SA6    YR
          SA7    YR+1          (YR) = SGN(A)*DSQRT(ONE+XR*XR)
*
*
          SA2    ONE           (X2) = +1.
*
          FX1    X2/X6         (X6/X7) = 1.D0 / YR
          FX4    X1*X6
          FX5    X2-X4
          DX4    X2-X4
          NX5    X5
          FX4    X4+X5
          DX5    X1*X6
          FX0    X1*X7
          FX4    X4-X5
          FX4    X4-X0
          FX0    X4/X6
          FX4    X0+X1
          DX5    X0+X1
          NX1    X4
          FX6    X1+X5
          DX7    X1+X5         (X6,X7) = (X2,) / (X6,X7)
          SA6    B3            DC=(X6,X7)
          SA7    B3-B7
*
          SA2    XR            (X2,X3) = XR
          SA3    XR+1
*
          FX4    X2*X7         (X6,X7) = DC * XR
          FX5    X3*X6
          FX4    X4+X5
          FX7    X2*X6
          DX5    X2*X6
          FX5    X4+X5
          FX6    X5+X7
          DX7    X5+X7         (X6,X7) = (X6,X7) * (X2,X3)
*
          SA6    B4            DS=(X6,X7)
          SA7    B4-B7
*
          SA2    B1
          SA3    B1-B7         (X2,X3)=DA
*
          SA4    YR            (X4,X5) = YR
          SA5    YR+1
*
          FX0    X3*X4         (X6,X7) = DA * YR
          FX1    X2*X5
          FX0    X0+X1
          FX3    X2*X4
          DX1    X2*X4
          FX1    X0+X1
          FX6    X1+X3
          DX7    X1+X3         (X6,X7) = (X2,X3) * (X4,X5)
*
          SA6    B1            DA=(X6,X7)
          SA7    B1-B7
*
          EQ     FIFTY         GO TO FIFTY
*
 THIRTY   MX6    0             (X6) = 0
          SA1    ONE           (X1) = +1.
          MX7    0             (X7) = 0
          SA6    B4            (DS) = (X6,X7)
          SA7    B4-B7         (DS) = (0.,0)
          BX6    X1            (X6) = X1
          SA7    B3-B7         DC = (X6,X7)
          SA6    B3
          EQ     FIFTY         GO TO FIFTY
*
 FORTY    MX6    0             (X6) = 0
          SA1    ONE           (X1) = +1.
          MX7    0             (X7) = 0
          SA6    B3            DC = (X6,X7)
          SA7    B3-B7
          BX6    X1            (X6) = +1.
          SA6    B4
          SA7    B4-B7         DS = (X6,X7)
*
          SA1    B2            (X1,X2) = DB
          SA2    B2-B7
          BX6    X1
          BX7    X2
          SA6    B1
          SA7    B1-B7         DA = (X1,X2)
*
 FIFTY    SA1    B3            (X1,) = DC
          ZR     X1,SEVENTY    IF(SNGL(DC) = 0)  GO TO SEVENTY
*
          BX6    X1            (X6) = X1
          AX1    59
          BX2    X6-X1         (X2,) = DABS(DC)
          SA3    B4            (X3,) = DS
          BX7    X3            (X7) = X3
          AX3    59
          BX4    X7-X3         (X4,) = DABS(DS)
*
          FX5    X4-X2         COMPARE UPPER HALVES:DABS(DC),DABS(DS)
          NX5    X5            MAKE SURE X5 DOES NOT CONTAIN A MINUS 0
          NG     X5,SIXTY      IF(DABS(DC) > ABS(DS)) GO TO SIXTY
*
          SA4    B3            (X4,X5) = DC
          SA5    B3-B7
          SA2    ONE           (X2) = +1.
          BX6    X4
          BX7    X5            (X6,X7) = DC
*
          FX1    X2/X6         (X6,X7) = 1.D0 / DC
          FX4    X1*X6
          FX5    X2-X4
          DX4    X2-X4
          NX5    X5
          FX4    X4+X5
          DX5    X1*X6
          FX0    X1*X7
          FX4    X4-X5
          FX4    X4-X0
          FX0    X4/X6
          FX4    X0+X1
          DX5    X0+X1
          NX1    X4
          FX6    X1+X5
          DX7    X1+X5         (X6,X7) = (X2,)/(X6,X7)
*
          SA6    B2            DB = 1.D0 / DC
          SA7    B2-B7
*
          EQ     OUT           GO TO OUT
*
 SIXTY    SA4    B4            (X4,X5) = DS
          SA5    B4-B7
          BX6    X4
          BX7    X5            (X6,X7) = (X4,X5)
          SA6    B2
          SA7    B2-B7         DB = (X6,X7)
*
          EQ     OUT           GO TO OUT
*
 SEVENTY  SA2    ONE           (X2) = +1.
          MX7    0             (X7) = 0.
          BX6    X2
          SA6    B2
          SA7    B2-B7         DB = (X6,X7)
*
 OUT      OUTFTN DROTG         RETURN
*
 ONE      DATA   17204000000000000000B
 XR       BSS    2
 XR2P1    BSS    2             TEMPORARY STORAGE FOR THE QUANTITY (1.+XR
 YR       BSS    2
 ADRDA    BSS    1
 ADRDB    BSS    1
 ADRDC    BSS    1
 ADRDS    BSS    1
*
          END
*DECK,SROT
          IDENT  SROT
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL SROT(N,SX,INCX,SY,INCY,SC,SS)
*
*         APPLY GIVENS REFLECTION MATRIX
*
*         APPLY 2X2 MATRIX  ( SC SS)  TO 2XN MATRIX  (SXI1  ... SXIN )
*                           (-SS SC)                 (SYI1  ... SYIN )
*
*         SXII  = SX(1 + (I-1)*INCX)  IF INCX .GE. 0
*               = SX(1 + (I-N)*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR SYII
*
*         SX( ),SY( )               SINGLE PRECISION
*         N,INCX,INCY               INTEGER TYPE
*         SC,SS                     SINGLE PRECISION
*
*         ROUNDED ARITHMETIC INSTRUCTIONS ARE USED
*
*         WRITTEN BY  RICHARD J. HANSON
*                     SANDIA LABORATORIES
*                     ALBUQUERQUE, NEW MEXICO
***       1 JUNE 77
*
          ENTRY  SROT
 SS       BSS    1
          VFD    42/4HSROT,18/7
*
 SROT     DATA   0
          INFTN  SROT,7      PROPER LINKAGE (RUN,FTN) MACRO.
          SA1    B1          (X1)=N
          SB7    1           (B7)=1
*
          SB1    X1          (B1)=N
          SB1    B1-B7       (B1)=N-1
*
          MI     B1,OUT      IF (N .LE. 0), QUIT
*
          SA5    SS          (X5)=LOC(SS)
          SA5    X5          (X5)=SS
*
          NZ     X5,APPLY    IF(SS.EQ.0..AND.SC.EQ.1.) QUIT.
          SA2    B6          (X2)=SC
          SA3    SONE        (X3)=1.
          RX2    X2-X3       (X2)=SC-1.
          NX2    X2          (X2)=NORM.(X2)
          ZR     X2,OUT      IF(SC.EQ.1.) QUIT.
 APPLY    SA1    B2          (X1)=SX(1)
          SA2    B3          (X2)=INCX
*
          SA3    B4          (X3)=SY(1)
          SA4    B5          (X4)=INCY
*
          ZR     B1,INCYNN   IF (N .EQ. 1) NO NEED TO TEST FOR NEG. INC.
          SX0    -B1         (X0)=-(N-1)
          SB2    X2          (B2)=INCX
          SB3    X4          (B3)=INCY
*
          GE     B2,INCXNN   IF (INCX .GE. 0) NO ADDRESS FIXUP NEEDED
          DX2    X0*X2       (X2)=-(N-1)*INCX
          SB7    A1          (B7)=LOC(SX(1))
          SA1    B7+X2       (X1)=SX(1+(1-N)*INCX),(A1)=LOC(X(1))
*
 INCXNN   GE     B3,INCYNN   IF (INCY .GE. 0) NO ADDRESS FIXUP NEEDED
          DX4    X0*X4       (X4)=-(N-1)*INCY
          SB7    A3          (B7)=LOC(SY(1))
          SA3    B7+X4       (X3)=SY(1+(1-N)*INCY),(A3)=LOC(Y(1))
*
 INCYNN   SA2    B6          (X2)=SC
          SB7    1           (B7)=1
          SB6    B7          (B6)=I=1
          BX0    X2          (X0)=SC
*
          SB1    B1-B7       (B1)=N-2
          GT     B6,B1,FIX   IF (I .GT. N-2) CLEAN-UP LOGIC
*
 LOOP     SA2    A1+B2       (X2)=SX(I+1)
          SA4    A3+B3       (X4)=SY(I+1)
          RX6    X3*X5       (X6)=SS*SY(I)
          RX7    X0*X3       (X7)=SC*SY(I)
          RX3    X0*X1       (X3)=SC*SX(I)
          RX1    X1*X5       (X1)=SS*SX(I)
*
          SB6    B6+B7       (B6)=I=I+1. INCREMENT I
          RX6    X3+X6       (X6)=SC*SX(I-1)+SS*SY(I-1)
          RX3    X4*X5       (X3)=SS*SY(I)
          RX7    X7-X1       (X7)=-SS*SX(I-1)+SC*SY(I-1)
          RX1    X0*X4       (X1)=SC*SY(I)
          RX4    X0*X2       (X4)=SC*SX(I)
          NX6    X6          (X6)=NORM.(X6)
          RX2    X2*X5       (X2)=SS*SX(I)
          NX7    X7          (X7)=NORM.(X7)
          NO     0           DEAD
          SA6    A1          SX(I-1)=(X6)
          NO     0           DEAD
          RX4    X3+X4       (X4)=SC*SX(I)+SS*SY(I)
          SA7    A3          SY(I-1)=(X7)
          SA3    A4+B3       (X3)=SY(I+1). NEXT ITERATION.
          RX2    X1-X2       (X2)=-SS*SX(I)+SC*SY(I)
          SA1    A2+B2       (X1)=SX(I+1). NEXT ITERATION.
          NX6    X4          (X6)=NORM.(X4)
          SB6    B6+B7       (B6)=I=I+1. INCREMENT I.
          NO     0           DEAD
          NX7    X2          (X7)=NORM(X2)
          SA6    A2          SX(I-1)=(X6)
          NO     2           DEAD
          SA7    A4          SY(I-1)=(X7)
          LE     B6,B1,LOOP  IF (I .LE. N-2) CONTINUE LOOP
 FIX      SB1    B1+B7       (B1)=N-1
          SB1    B1+B7       (B1)=N
 CL       RX6    X3*X5       (X6)=SS*SY(I)
          RX7    X0*X3       (X7)=SC*SY(I)
          RX3    X0*X1       (X3)=SC*SX(I)
          RX1    X1*X5       (X1)=SS*SX(I)
*
          SB6    B6+B7       (B6)=I=I+1. INCREMENT I.
          RX6    X3+X6       (X6)=SC*SX(I-1)+SS*SY(I-1)
          RX7    X7-X1       (X7)=-SS*SX(I-1)+SC*SY(I-1)
*
          NX6    X6          (X6)=NORM.(X6)
          NX7    X7          (X7)=NORM.(X7)
*
          SA6    A1          SX(I-1)=(X6)
          SA7    A3          SY(I-1)=(X7)
*
          GT     B6,B1,OUT   IF (I .GT. N), QUIT
          SA3    A3+B3       (X3)=SY(I)
          SA1    A1+B2       (X1)=SX(I)
          JP     CL          ONE COMP. REMAINS.
 OUT      OUTFTN SROT
 SONE     DATA   1.0
*         END    SROT
          END
*DECK,DROT
          IDENT  DROT
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL DROT(N,DX,INCX,DY,INCY,DC,DS)
*
*         APPLY GIVENS REFLECTION MATRIX
*
*         APPLY 2X2 MATRIX  ( DC DS)  TO 2XN MATRIX  (DXI1  ... DXIN )
*                           (-DS DC)                 (DYI1  ... DYIN )
*
*         DXII  = DX(1 + (I-N)*2*INCX)  IF INCX .GE. 0
*               = DX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR DYII
*
*         DX( ),DY( )               DOUBLE PRECISION
*         N,INCX,INCY               INTEGER TYPE
*         DC,DS                     DOUBLE PRECISION
*
*         WRITTEN BY  DAVID R. KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  DROT
 ARG7     BSS    1
          VFD    42/4HDROT,18/7
*
 DROT     DATA   0             ENTRY/EXIT
          INFTN  DROT,7
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          SB1    X1+B7         (B1) = N-1
*
          SA2    B3            (X2) = INCX
          NG     B1,OUT        IF N .LE. 0 , GO TO OUT
          SA3    ARG7          (X3) = LOC(DS)
          SA3    X3            (X3,) = DS
          NZ     X3,DROT5      IF DS.NE.0.0EE0, GO TO DROT5
*
          SA3    B6            (X3,X4) = DC
          SA4    B6-B7
          SA1    DONE          (X1,X5) = 1.0EE0
          SA5    A1-B7
*
          FX4    X4-X5
          FX5    X3-X1
          DX3    X3-X1
          NX5    X5
          FX3    X3+X4
          NX3    X3
          FX5    X3+X5
          ZR     X5,OUT        IF DC.EQ.1.0EE0.AND.DS.EQ.0.0EE0, GOTO OU
*
 DROT5    SA3    B5            (X3) = INCX
          SX1    -B1           (X1) = -(N-1)
          LX2    1             INCX = 2*INCX
          IX3    X3+X3         INCY = 2*INCY
          SB3    X2            (B3) = INCX
          SB5    X3            (B5) = INCY
*
          GT     B3,DROT10     IF INCX .GT. 0 , GO TO DROT10
          ZR     B3,OUT        IF INCX .EQ. 0 , GO TO OUT
          DX0    X1*X2         LOC(DXI1 ) = LOC(DX) - (N-1)*INCX
          SB2    X0+B2         (B2) = LOC(DXI1 )
*
 DROT10   GT     B5,DROT20     IF INCY .GT. 0, GO TO DROT20
          ZR     B5,OUT        IF INCY .EQ. 0 , GO TO OUT
          DX0    X1*X3         LOC(DYI1 ) = LOC(DY) - (N-1)*INCY
          SB4    X0+B4         (B4) = LOC(DYI1 )
*
 DROT20   SA5    ARG7
          SA0    X5            (A0) = LOC(DS)
          SB1    B1-B7         (B1) = N
*
 LOOP     SA1    A0            (X1,X2) = DS
          SA2    A0-B7
*
          SA3    B4            (X3,X4) = DYII
          SA4    B4-B7
*
          FX5    X2*X3         (X6,X7) = DS*DYII
          FX0    X1*X4
          FX5    X0+X5
          FX4    X1*X3
          DX0    X1*X3
          FX5    X0+X5
          FX6    X4+X5
          DX7    X4+X5
*
          SA1    B2            (X1,X2) = DXII
          SA2    B2-B7
*
          SA3    B6            (X3,X4) = DC
          SA4    B6-B7
*
          FX5    X2*X3         (X0,X3) = DC*DXII
          FX0    X1*X4
          FX5    X0+X5
          FX4    X1*X3
          DX0    X1*X3
          FX5    X0+X5
          FX0    X4+X5
          DX3    X4+X5
*
          FX4    X6+X0         (X6,X7) = (X6,X7)+(X0,X3)
          DX5    X6+X0
          FX0    X7+X3
          NX4    X4
          FX3    X0+X5
          FX0    X3+X4
          NX5    X0
          DX3    X3+X4
          NX4    X3
          FX6    X4+X5
          DX7    X4+X5
*
          SA6    DW              DW = (X6,X7)
          SA7    DW+1
*
          SA3    A0            (X3,X4) = DS
          SA4    A0-B7
*
          FX5    X2*X3         (X6,X7) = DS*DXII
          FX0    X1*X4
          FX5    X0+X5
          FX4    X1*X3
          DX0    X1*X3
          FX5    X0+X5
          FX6    X4+X5
          DX7    X4+X5
*
          SA1    B6            (X1,X2) = DC
          SA2    B6-B7
          SA3    B4            (X3,X4) = DYII
          SA4    B4-B7
*
          FX5    X2*X3         (X0,X2) = DC*DYII
          FX0    X1*X4
          FX5    X0+X5
          FX4    X1*X3
          DX0    X1*X3
          FX5    X0+X5
          FX0    X4+X5
          DX2    X4+X5
*
          FX4    X0-X6         (X6,X7) = (X0,X2)-(X6,X7)
          DX5    X0-X6
          FX0    X2-X7
          NX4    X4
          FX2    X0+X5
          FX0    X2+X4
          NX5    X0
          DX2    X2+X4
          NX4    X2
          FX6    X4+X5
          DX7    X4+X5
*
          SA6    B4            DYII  = (X6,X7)
          SA7    B4-B7
*
          SB1    B1+B7         COUNT TERM
          SA1    DW
          SA2    DW+1
          BX6    X1
          BX7    X2
          SA6    B2
          SA7    B2-B7         DXII  = DW
*
          SB2    B2+B3         (B2) = LOC(DXII+1 )
          SB4    B4+B5         (B4) = LOC(DYII+1 )
*
          NZ     B1,LOOP       IF I .NE. N, LOOP
*
 OUT      OUTFTN DROT          RETURN
*
 DONE     DATA   1.0EE0
 DW       BSS    2
*
          END
*DECK,SROTMG
          IDENT  SROTMG
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL SROTMG(SD1,SD2,SB1,SB2,SPARAM)
*
*         CONSTRUCT THE TWO-MULTIPLY,TWO-ADD,NO-SQUARE-RO0T
*         GIVENS ROTATION
*
*
*         THIS SUBROUTINE STORES VALUES IN SPARAM( )
*         DEFINING THE MATRIX H
*
*         SPARAM(1) = FLAG , INDICATES THE FORM OF THE MATRIX H
*         SPARAM(2) = H11
*         SPARAM(3) = H21
*         SPARAM(4) = H12
*         SPARAM(5) = H22
*
*         THE FLAG VALUES AND THE CORRESPONDING FORMS OF THE MATRIX H
*         -2. (1 0)   -1. (H11 H12)   0. ( 1 H12)   1. (H11 1 )
*             (0 1)       (H21 H22)      (H21 1 )      (-1 H22)
*
*         SD1,SD2,SB1,SB2           SINGLE PRECISION
*         SPARAM( )                 SINGLE PRECISION
*
*         THIS ALGORITHM ASSUMES THAT THE INPUT VALUE OF SD1 IS
*         POSITIVE OR ZERO BUT NON-NEGATIVE. THE VALUE OF SD2 IS
*         UNRESTRICTED.
*
*         ROUNDED ARITHMETIC INSTRUCTIONS ARE USED
*
*         WRITTEN BY  DAVID R. KINCAID AND JAMES SULLIVAN
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  SROTMG
          VFD    42/6HSROTMG,18/5
*
 SROTMG   DATA   0             ENTRY/EXIT
          INFTN  SROTMG,5
          SA1    B1            (X1) = SD1
          SA3    B3            (X3) = SB1
          RX0    X1*X3         (X0) = P1 = SD1*SB1
          SA2    B2            (X2) = SD2
          SA4    B4            (X4) = SB2
          RX5    X2*X4         (X5) = P2 = SD2*SB2
          RX6    X0*X3         (X6) = P1*SB1
          RX7    X5*X4         (X7) = P2*SB2
*
          BX1    X6
          AX6    59
          BX6    X6-X1         (X6) = ABS(P1*SB1)
*
          BX2    X7
          AX7    59
          BX7    X7-X2         (X7) = ABS(P2*SB2)
*
          RX6    X7-X6
          NX6    X6
          NG     X6,TWELVE     IF( ABS(P1*SB1) .GT. ABS(P2*SB2) )
*                                     GO TO 12
*
          ZR     X2,FOUR
          NG     X2,SIXTN      IF( P2*SB2 ) 16,4,10
*
*
          RX7    X3/X4         (X7) = SB1/SB2      ITEN
          SA1    B1            (X1) = SD1
          SA2    B2            (X2) = SD2
          RX6    X0/X5         (X6) = P1/P2
          SA7    B5+4          SPARAM(5) = (X7)
          SA6    B5+1          SPARAM(2) = (X6)
          RX0    X6*X7         (X0) = SPARAM(2)*SPARAM(5)
          SA5    UNIT          (X5) = 1.0
          RX0    X5+X0         (X0) = 1.0 + SPARAM(2)*SPARAM(5) = U
          NX0    X0
          BX7    X5            (X7) = X5
          RX5    X5/X0         (X5) = 1./U
          SA7    B5            SPARAM(1) = 1.0
          RX7    X4*X0         (X7) = SB2*(X0)
          SA7    B3            SB1 = (X7)
          BX3    X7            (X3) = SB1
          RX6    X2*X5         (X6) = SD2*(X5)
          RX7    X1*X5         (X7) = SD1*(X5)
          SA6    B1            SD1 = (X6)
          SA7    B2            SD2 = (X7)
          BX1    X6            (X1) = SD1
          BX2    X7            (X2) = SD2
          EQ     TWENTY4       GO TO 24
*
FOUR      SA5    RTWO
          BX6    -X5
          SA6    B5            SPARAM(1) = -2.0
*
          EQ     OUT           GO TO OUT
*
*
 TWELVE   RX7    X4/X3         (X7) = SB2/SB1
          SA1    B1            (X1) = SD1
          SA2    B2            (X2) = SD2
          RX6    X5/X0         (X6) = P2/P1
          BX7    -X7           (X7) = -SB2/SB1
          SA7    B5+2          SPARAM(3) = (X7)
          SA6    B5+3          SPARAM(4) = (X6)
          RX0    X6*X7         (X0) = SPARAM(4)*SPARAM(3)
          SA5    UNIT          (X5) = 1.0
          RX5    X5-X0         (X0) = 1.0 - SPARAM(4)*SPARAM(3) = U
          NX0    X5
*
          SA5    TOL           (X5) = TOL
          RX5    X5-X0
          NX5    X5
          PL     X5,SIXTN     IF( U .LE. TOL ) GO TO 16
*
*                              HERE WHEN U IS ZERO OR NEARLY ZERO.
*                              ALSO WHEN SD1 IS NEGATIVE AND
*                              ABS(SD1*SB1**1) .LE. ABS(SD2*SB2**2)
*                              SINCE IN SUCH A CASE U SHOULD BE SMALL.
*
          SA5    UNIT          (X5) = 1.0
          RX5    X5/X0         (X5) = 1./U
          RX7    X3*X0         (X7) = SB1*U
          SA7    B3            SB1 = (X7)
          MX6    0             (X6) = 0.0
          SA6    B5            SPARAM(1) = 0.0
          BX3    X7            (X3) = SB1
          RX6    X1*X5         (X6) = SD1*(X5)
          RX7    X2*X5         (X7) = SD2*(X5)
          SA6    A1            SD1 = (X6)
          SA7    A2            SD2 = (X7)
          BX1    X6            (X1) = SD1
          BX2    X7            (X2) = SD2
*
          EQ     TWENTY4       RETURN
*
 SIXTN    MX7    0             (X7) = 0.0
          SA5    UNIT          (X5) = -1.0
          BX6    -X5
          SA6    B5            SPARAM(1) = -1.0
          SA7    B5+1          SPARAM(2) = 0.0
          MX6    0             (X6) = 0.0
          SA6    B5+2          SPARAM(3) = 0.0
          SA7    B5+3          SPARAM(4) = 0.0
          SA6    B5+4          SPARAM(5) = 0.0
          SA7    B1            SD1 = 0.0
          SA6    B2            SD2 = 0.0
          SA7    B3            SB1 = 0.0
*
          EQ     OUT           GO TO OUT
*
 TWENTY4  BX6    X1
          SA5    BIGINV
          AX6    59
          BX6    X6-X1
          RX5    X5-X6
          NX5    X5
          NG     X5,THIRTY6
          ZR     X1,FOURTY8
          SA5    B5
          ZR     X5,A84
          NG     X5,A32
          SA5    UNIT
          BX6    X5
          BX7    -X5
          SA6    B5+3
          SA7    B5+2
          EQ     A92
 A84      SA5    UNIT
          BX6    X5
          BX7    X5
          SA6    B5+1
          SA7    B5+4
          BX7    -X5
 A92      SA7    B5
 A32      SA5    SQRBIG2
          RX6    X1*X5
          SA5    SQRBIGI
          BX1    X6
          SA6    B1
          SA4    B5+1
          RX6    X3*X5
          RX7    X4*X5
          SA6    B3
          SA7    B5+1
          BX3    X6
          SA4    B5+3
          RX6    X4*X5
          SA6    B5+3
          EQ     TWENTY4
 THIRTY6  BX6    X1
          SA5    BIG
          AX6    59
          BX6    X6-X1
          RX5    X6-X5
          NX5    X5
          NG     X5,FOURTY8
          SA5    B5
          ZR     X5,B84
          NG     X5,B32
          SA5    UNIT
          BX6    X5
          BX7    -X5
          SA6    B5+3
          SA7    B5+2
          EQ     B92
 B84      SA5    UNIT
          BX6    X5
          BX7    X5
          SA6    B5+1
          SA7    B5+4
          BX7    -X5
 B92      SA7    B5
 B32      SA5    SQRBI2I
          RX6    X1*X5
          SA5    SQRBIG
          BX1    X6
          SA6    B1
          SA4    B5+1
          RX6    X3*X5
          RX7    X4*X5
          SA6    B3
          SA7    B5+1
          BX3    X6
          SA4    B5+3
          RX6    X4*X5
          SA6    B5+3
          EQ     THIRTY6
 FOURTY8  BX4    X2
          SA5    BIGINV
          AX4    59
          BX4    X4-X2
          RX5    X5-X4
          NX5    X5
          NG     X5,SIXTY
          ZR     X2,OUT
          SA5    B5
          ZR     X5,C84
          NG     X5,C32
          SA5    UNIT
          BX6    X5
          BX7    -X5
          SA6    B5+3
          SA7    B5+2
          EQ     C92
 C84      SA5    UNIT
          BX6    X5
          BX7    X5
          SA6    B5+1
          SA7    B5+4
          BX7    -X5
 C92      SA7    B5
 C32      SA5    SQRBIG2
          RX6    X2*X5
          SA5    SQRBIGI
          BX2    X6
          SA6    B2
          SA4    B5+2
          RX7    X4*X5
          SA7    B5+2
          SA4    B5+4
          RX6    X4*X5
          SA6    B5+4
          EQ     FOURTY8
 SIXTY    BX4    X2
          SA5    BIG
          AX4    59
          BX4    X4-X2
          RX5    X4-X5
          NX5    X5
          NG     X5,OUT
          SA5    B5
          ZR     X5,D84
          NG     X5,D32
          SA5    UNIT
          BX6    X5
          BX7    -X5
          SA6    B5+3
          SA7    B5+2
          EQ     D92
 D84      SA5    UNIT
          BX6    X5
          BX7    X5
          SA6    B5+1
          SA7    B5+4
          BX7    -X5
 D92      SA7    B5
 D32      SA5    SQRBI2I
          RX6    X2*X5
          SA5    SQRBIG
          BX2    X6
          SA6    B2
          SA4    B5+2
          RX7    X4*X5
          SA7    B5+2
          SA4    B5+4
          RX6    X4*X5
          SA6    B5+4
          EQ     SIXTY
 OUT      OUTFTN SROTMG        RETURN
*
 BIG      DATA   1.67772E7
 BIGINV   DATA   5.96046E-8
 RTWO     DATA   2.0
 SQRBIG   DATA   4096.0
 SQRBIGI  DATA   17044000000000000000B
 SQRBIG2  DATA   17504000000000000000B
 SQRBI2I  DATA   16704000000000000000B
 TOL      DATA   0.0
 UNIT     DATA   1.0
*
          END
*DECK,DROTMG
          IDENT  DROTMG
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL DROTMG(DD1,DD2,DB1,DB2,DPARAM)
*
*         CONSTRUCT THE TWO-MULTIPLY,TWO-ADD,NO-SQUARE-RO0T
*         GIVENS ROTATION
*
*
*         THIS SUBROUTINE STORES VALUES IN DPARAM( )
*         DEFINING THE MATRIX H
*
*         DPARAM(1) = FLAG , INDICATES THE FORM OF THE MATRIX H
*         DPARAM(2) = H11
*         DPARAM(3) = H21
*         DPARAM(4) = H12
*         DPARAM(5) = H22
*
*         THE FLAG VALUES AND THE CORRESPONDING FORMS OF THE MATRIX H
*         -2. (1 0)   -1. (H11 H12)   0. ( 1 H12)   1. (H11 1 )
*             (0 1)       (H21 H22)      (H21 1 )      (-1 H22)
*
*         DD1,DD2,DB1,DB2           DOUBLE PRECISION
*         DPARAM( )                 DOUBLE PRECISION
*
*         THIS ALGORITHM ASSUMES THAT THE INPUT VALUE OF DD1 IS
*         POSITIVE OR ZERO BUT NON-NEGATIVE. THE VALUE OF DD2 IS
*         UNRESTRICTED.
*
*         WRITTEN BY  DAVID R. KINCAID AND JAMES SULLIVAN
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  DROTMG
          VFD    42/6HDROTMG,18/5
*
 DROTMG   DATA   0             ENTRY/EXIT
          INFTN  DROTMG,5
          SA1    B1            (X1,X2) = DD1
          SA2    B1+1
          SA3    B3            (X3,X4) = DB1
          SA4    B3+1
*
          FX7    X2*X3         (X6,X7) = DD1 * DB1
          FX6    X1*X4
          FX7    X6+X7
          DX6    X1*X3
          FX0    X1*X3
          FX7    X6+X7
          FX6    X0+X7
          DX7    X0+X7         (X6,X7) = (X1,X2) * (X3,X4)
*
          SA6    P1           (P1) = (X6,X7)
          SA7    P1+1         (P1) = DD1 * DB1
*
          FX1    X4*X6         (X0,X1) = P1 * DB1
          FX2    X3*X7
          FX1    X1+X2
          DX0    X3*X6
          FX2    X3*X6
          FX1    X0+X1
          FX0    X1+X2
          DX1    X1+X2         (X0,X1) = (X3,X4) * (X6,X7)
*
          BX2    X0
          AX2    59
          BX0    X2-X0
          BX1    X2-X1         (X0,X1) = DABS( P1*DB1 )
*
          SA2    B2            (X2,X3) = DD2
          SA3    B2+1
          SA4    B4            (X4,X5) = DB2
          SA5    B4+1
*
          FX7    X3*X4         (X6,X7) = DD2 * DB2
          FX6    X2*X5
          FX7    X6+X7
          FX3    X2*X4
          DX6    X2*X4
          FX7    X6+X7
          FX6    X3+X7
          DX7    X3+X7         (X6,X7) = (X2,X3) * (X4,X5)
*
          SA6    P2           (P2) = (X6,X7)
          SA7    P2+1         (P2) = DD2 * DB2
*
          FX2    X5*X6         (X2,X3) = P2 * DB2
          FX3    X4*X7
          FX2    X2+X3
          FX7    X4*X6
          DX3    X4*X6
          FX3    X2+X3
          FX2    X3+X7
          DX3    X3+X7         (X2,X3) = (X4,X5) * (X6,X7)
*
          BX6    X2
          BX7    X3
          SA6    TEMP
          SA7    TEMP+1        TEMP = P2*DB2
*
          AX6    59
          BX2    X6-X2
          BX3    X6-X3         (X2,X3) = DABS( P2*DB2 )
*
          FX6    X2-X0         COMPUTE DABS(P2*DB2) - DABS(P1*DB1).
          DX7    X2-X0
          FX2    X3-X1
          NX6    X6
          FX0    X2+X7
          FX2    X0+X6
          NX7    X2
          DX0    X0+X6
          NX6    X0
          FX2    X6+X7
          DX3    X6+X7         (X2,X3) = (X2,X3) - (X0,X1)
*
*
          NG     X2,TWELVE     IF( DABS(P1*DB1) .GT. DABS(P2*DB2) )
*                                    GO TO TWELVE
*
          SA2    TEMP          (X2,X3) = P2*DB2
          SA3    TEMP+1
          ZR     X2,FOUR
          NG     X2,SIXTN      IF( P2*DB2 ) SIXTN,FOUR,TEN
*
          SA2    B3            (X2,X3) = DB1      ITEN
          SA3    B3+1
          SA4    B4            (X4,X5) = DB2
          SA5    B4+1
*
          FX1    X2/X4         (X6,X7) = DB1 / DB2
          FX6    X1*X4
          FX7    X2-X6
          DX6    X2-X6
          NX7    X7
          FX6    X6+X7
          DX7    X1*X4
          FX0    X1*X5
          FX6    X3+X6
          FX6    X6-X7
          FX6    X6-X0
          FX0    X6/X4
          FX6    X0+X1
          DX7    X0+X1
          NX1    X6
          FX6    X1+X7
          DX7    X1+X7         (X6,X7) = (X2,X3) / (X4,X5)
*
          SA6    B5+8           (DPARAM(5)) = (X6,X7)
          SA7    B5+9           (DPARAM(5)) = DB1 / DB2
*
          SA2    P1           (X2,X3) = P1
          SA3    P1+1
          SA4    P2           (X4,X5) = P2
          SA5    P2+1
*
          FX1    X2/X4         (X6,X7) = P1 / P2
          FX6    X1*X4
          FX7    X2-X6
          DX6    X2-X6
          NX7    X7
          FX6    X6+X7
          DX7    X1*X4
          FX0    X1*X5
          FX6    X3+X6
          FX6    X6-X7
          FX6    X6-X0
          FX0    X6/X4
          FX6    X0+X1
          DX7    X0+X1
          NX1    X6
          FX6    X1+X7
          DX7    X1+X7         (X6,X7) = (X2,X3) / (X4,X5)
*
          SA6    B5+2         (DPARAM(2)) = (X6,X7)
          SA7    B5+3         (DPARAM(2)) = P1 / P2
*
          SA4    B5+8                (X4,X5) = B5+8
          SA5    B5+9
*
          FX1    X4*X7         (X1,X2) = DPARAM(2) * DPARAM(5)
          FX2    X5*X6
          FX1    X1+X2
          DX2    X4*X6
          FX0    X4*X6
          FX1    X1+X2
          DX2    X0+X1
          FX1    X0+X1         (X1,X2) = (X4,X5) * (X6,X7)
*
          SA3    UNIT          (X3) = +1.
*
          FX6    X1+X3         (X4,X5) = 1.D0 + (DPARAM(2)*DPARAM(5))
          DX7    X1+X3
          NX6    X6
          FX5    X2+X7
          FX4    X5+X6
          NX7    X4
          DX5    X5+X6
          NX6    X5
          FX4    X6+X7
          DX5    X6+X7         (X4,X5) = (X3,0) + (X1,X2)    IU
*
          SA1    B4            (X1,X2) = DB2
          SA2    B4+1
*
          FX6    X1*X5         (X6,X7) = DB2 * U
          FX7    X2*X4
          FX6    X6+X7
          DX7    X1*X4
          FX0    X1*X4
          FX7    X6+X7
          FX6    X0+X7
          DX7    X0+X7         (X6,X7) = (X1,X2) * (X4,X5)
*
          SA6    B3            (DB1) = (X6,X7)
          SA7    B3+1          (DB1) = DB2 * U
*
          FX7    X3/X4         (X0,X1) = 1.D0 / U
          FX0    X4*X7
          FX1    X3-X0
          DX0    X3-X0
          NX1    X1
          FX0    X0+X1
          DX1    X4*X7
          FX6    X5*X7
          FX0    X0-X1
          FX0    X0-X6
          FX6    X0/X4
          FX0    X6+X7
          DX1    X6+X7
          NX7    X0
          FX0    X1+X7
          DX1    X1+X7         (X0,X1) = (X3,0) / (X4,X5)
*
          SA2    B2            (X2,X3) = DD2
          SA3    B2+1
*
          FX6    X1*X2         (X6,X7) = DD2 * (1.D0/U)
          FX7    X0*X3
          FX6    X6+X7
          DX7    X0*X2
          FX4    X0*X2
          FX7    X6+X7
          FX6    X4+X7
          DX7    X4+X7         (X6,X7) = (X0,X1) * (X2,X3)   IZ
*
          SA4    B1            (X4,X5) = DD1
          SA5    B1+1
*
          SA6    A4            (DD1) = (X6,X7)
          SA7    A5            (DD1) = DD2 * (1.D0/U) = Z
*
          FX2    X1*X4         (X6,X7) = DD1 * (1.D0/U)
          FX3    X0*X5
          FX2    X2+X3
          DX3    X0*X4
          FX7    X0*X4
          FX3    X2+X3
          FX6    X3+X7
          DX7    X3+X7         (X6,X7) = (X0,X1) * (X4,X5)
*
          SA6    A2            (DD2) = (X6,X7)
          SA7    A3            (DD2) = DD1 * (1.D0/U)
          SA1    UNIT
          MX7    0
          BX6    X1
          SA7    B5+1
          SA6    B5
*
          EQ     TWENTY4       GO TO TWENTY4
 FOUR     SA1    RTWO          (X1) = 2.0
          MX7    0             (X7) = 0.0
          BX6    -X1           (X6) = -(X1)
          SA7    B5+1
          SA6    B5            DPARAM(1) = -2.0
          EQ     OUT           GO TO OUT
*
*
 TWELVE   SA2    B3            (X2,X3) = DB1
          SA3    B3+1
          SA4    B4            (X4,X5) = DB2
          SA5    B4+1
*
          FX1    X4/X2         (X6,X7) = DB2 / DB1
          FX6    X1*X2
          FX7    X4-X6
          DX6    X4-X6
          NX7    X7
          FX6    X6+X7
          DX7    X1*X2
          FX0    X1*X3
          FX6    X5+X6
          FX6    X6-X7
          FX6    X6-X0
          FX0    X6/X2
          FX6    X0+X1
          DX7    X0+X1
          NX1    X6
          FX6    X1+X7
          DX7    X1+X7         (X6,X7) = (X4,X5) / (X2,X3)
*
          BX6    -X6           (X6,X7) = -DB2/DB1
          BX7    -X7
*
          SA6    B5+4     (DPARAM(3)) = (X6,X7)
          SA7    B5+5                 = - DB2 / DB1
*
          SA2    P2           (X2,X3) = P2
          SA3    P2+1
          SA4    P1           (X4,X5) = P1
          SA5    P1+1
*
          FX1    X2/X4         (X6,X7) = P2 / P1
          FX6    X1*X4
          FX7    X2-X6
          DX6    X2-X6
          NX7    X7
          FX6    X6+X7
          DX7    X1*X4
          FX0    X1*X5
          FX6    X3+X6
          FX6    X6-X7
          FX6    X6-X0
          FX0    X6/X4
          FX6    X0+X1
          DX7    X0+X1
          NX1    X6
          FX6    X1+X7
          DX7    X1+X7         (X6,X7) = (X2,X3) / (X4,X5)
*
          SA6    B5+6     (DPARAM(4) = (X6,X7)
          SA7    B5+7     (DPARAM(4) = P2 / P1
*
          SA4    B5+4     (X4,X5) = DPARAM(3)
          SA5    B5+5
*
          FX1    X4*X7         (X1,X2) = DPARAM(4) * DPARAM(3)
          FX2    X5*X6
          FX1    X1+X2
          DX2    X4*X6
          FX0    X4*X6
          FX1    X1+X2
          DX2    X0+X1
          FX1    X0+X1         (X1,X2) = (X4,X5) * (X6,X7)
*
          SA3    UNIT          (X3) = +1.
*
          FX6    X3-X1         (X4,X5) = 1.D0 - (DPARAM(4)*DPARAM(3))
          DX7    X3-X1
          NX6    X6
          FX5    X7-X2
          FX4    X5+X6
          NX7    X4
          DX5    X5+X6
          NX6    X5
          FX4    X6+X7
          DX5    X6+X7         (X4,X5) = (X3,0) - (X1,X2)    IU
*
* INSERT IF(U .LE. TOL) GO TO 16  HERE
          ZR     X4,SIXTN
          SA1    B3            (X1,X2) = DB1
          SA2    B3+1
*
          FX6    X1*X5         (X6,X7) = DB1 * U
          FX7    X2*X4
          FX6    X6+X7
          DX7    X1*X4
          FX0    X1*X4
          FX7    X6+X7
          FX6    X0+X7
          DX7    X0+X7         (X6,X7) = (X1,X2) * (X4,X5)
*
          SA6    A1            (DB1) = (X6,X7)
          SA7    A2            (DB1) = DB1 * U
*
          FX7    X3/X4         (X0,X1) = 1.D0 / U
          FX0    X4*X7
          FX1    X3-X0
          DX0    X3-X0
          NX1    X1
          FX0    X0+X1
          DX1    X4*X7
          FX6    X5*X7
          FX0    X0-X1
          FX0    X0-X6
          FX6    X0/X4
          FX0    X6+X7
          DX1    X6+X7
          NX7    X0
          FX0    X1+X7
          DX1    X1+X7         (X0,X1) = (X3,0) / (X4,X5)
*
          SA2    B1            (X2,X3) = DD1
          SA3    B1+1
*
          FX6    X1*X2         (X6,X7) = DD1 * (1.D0/U)
          FX7    X0*X3
          FX6    X6+X7
          DX7    X0*X2
          FX4    X0*X2
          FX7    X6+X7
          FX6    X4+X7
          DX7    X4+X7         (X6,X7) = (X0,X1) * (X2,X3)
*
          SA6    A2            (DD1) = (X6,X7)
          SA7    A3            (DD1) = DD1 / U
*
          SA4    B2            (X4,X5) = DD2
          SA5    B2+1
*
          FX6    X1*X4         (X6,X7) = DD2 * (1.D0/U)
          FX7    X0*X5
          FX6    X6+X7
          DX7    X0*X4
          FX2    X0*X4
          FX7    X6+X7
          FX6    X2+X7
          DX7    X2+X7         (X6,X7) = (X0,X1) * (X4,X5)
*
          SA6    A4            (DD2) = (X6,X7)
          SA7    A5            (DD2) = DD2 / U
          MX6    0
          BX7    X6
          SA6    B5
          SA7    B5+1
*
          EQ     TWENTY4
*
*
*
*                              HERE WHEN U IS ZERO OR NEARLY ZERO.
*                              ALSO WHEN D1 IS NEGATIVE AND
*                              DABS(D1*B1**2) .LE. DABS(D2*B2**2)
*                              SINCE IN SUCH A CASE U SHOULD BE SMALL.
*
*
 SIXTN    SA5    UNIT          (X5) = +1.0
          MX4    0             (X4) = 0.0
          BX7    X4            (X7) = 0.0
          BX6    -X5           (X6) = -X5
          SA6    B5            (SPARAM(1)) = (X6,X7) = -1.0D
          SA7    B5+1
          BX6    X7            (X6,X7) = 0.0D
          SA6    B5+2          (SPARAM(2)) = 0.0D
          SA7    B5+3
          SA6    B5+4          (SPARAM(3)) = 0.0D
          SA7    B5+5
          SA6    B5+6          (SPARAM(4)) = 0.0D
          SA7    B5+7
          SA6    B5+8          (SPARAM(5)) = 0.0D
          SA7    B5+9
*
          SA6    B1            DD1 = 0.0D
          SA7    B1+1
          SA6    B2            DD2 = 0.0D
          SA7    B2+1
          SA6    B3            DB1 = 0.0D
          SA7    B3+1
*
          EQ     OUT           GO TO OUT
*
*
*                              HERE TO RESCALE IF NECESSARY TO KEEP
*                              DD1 AND DD2 BETWEEN  BIG AND 1/BIG
*                              IF NONZERO
*
*
 TWENTY4  SA3    B1            (X3,X4) = DD1
          SA4    B1+1
          SA5    BIGINV        (X5,) = BIGINV
          BX0    X3
          AX0    59
          BX0    X0-X3         (X0,) = DABS(DD1)
          FX0    X5-X0         IF ( DABS(DD1) .GT. BIGINV ) GO TO 36
          NX0    X0
          NG     X0,THIRTY6
          ZR     X3,FOURTY8    IF (DD1) 28,48,28
          SA1    B5            (X1,) = DPARAM(1)          I28
          SA2    UNIT
          ZR     X1,A84        IF (DPARAM(1)) 96,84,88(A)
          NG     X1,A96
          BX6    X2                                        IA88
          MX7    0
          SA6    B5+6          DPARAM(4) = 1.0
          SA7    B5+7
          BX6    -X6
          SA6    B5+4          DPARAM(3) = -1.0
          SA7    B5+5
          EQ     A92           GO TO 92(A)
 A84      BX6    X2
          MX7    0
          SA6    B5+2          DPARAM(2) = 1.0
          SA7    B5+3
          SA6    B5+8          DPARAM(5) = 1.0
          SA7    B5+9
          BX6    -X6
 A92      SA7    B5+1          DPARAM(1) = -1.0
          SA6    B5
 A96      SA5    BIG           (X5,) = BIG
          FX1    X4*X5         DD1 = DD1 * (SQRBIG*SQRBIG)
          DX7    X3*X5
          FX6    X3*X5
          FX1    X1+X7
          DX7    X1+X6
          FX6    X1+X6         (X6,X7) = (X3,X4) * (X5,0)
          SA7    B1+1          DD1 = (X6,X7)
          SA6    B1
          SA2    B3            (X2,X3) = DB1
          SA3    B3+1
          SA4    SQRBIGI       (X4,) = SQRBIGI
          FX1    X3*X4         (X6,X7) = DB1/SQRBIG
          DX7    X2*X4
          FX6    X2*X4
          FX1    X1+X7
          DX7    X1+X6
          FX6    X1+X6         (X6,X7) = (X2,X3) * (X4,0)
          SA7    B3+1          DB1 = (X6,X7)
          SA6    B3
          SA2    B5+2          (X2,X3) = DPARAM(2)
          SA3    B5+3
          FX1    X3*X4         (X6,X7) = DPARAM(2)/SQRBIG
          DX7    X2*X4
          FX6    X2*X4
          FX1    X1+X7
          DX7    X1+X6
          FX6    X1+X6         (X6,X7) = (X2,X3) * (X4,0)
          SA7    B5+3          DPARAM(2) = (X6,X7)
          SA6    B5+2
          SA2    B5+6          (X2,X3) = DPARAM(4)
          SA3    B5+7
          FX1    X3*X4         (X6,X7) = DPARAM(4)/SQRBIG
          DX7    X2*X4
          FX6    X2*X4
          FX1    X1+X7
          DX7    X1+X6
          FX6    X1+X6         (X6,X7) = (X2,X3) * (X4,0)
          SA7    B5+7          DPARAM(4) = (X6,X7)
          SA6    B5+6
          EQ     TWENTY4       GO TO 24
 THIRTY6  SA3    B1            (X3,X4) = DD1
          SA4    B1+1
          SA5    BIG           (X5,) = BIG
          BX0    X3
          AX0    59
          BX0    X0-X3         (X0,) = DABS(DD1)
          FX0    X0-X5         IF ( DABS(DD1) .LT. BIG ) GO TO 48
          NX0    X0
          NG     X0,FOURTY8
          SA1    B5            (X1,) = DPARAM(1)
          SA2    UNIT
          ZR     X1,B84        IF (DPARAM(1)) 96,84,88(B)
          NG     X1,B96
          BX6    X2                                        IB88
          MX7    0
          SA6    B5+6          DPARAM(4) = 1.0
          SA7    B5+7
          BX6    -X6
          SA6    B5+4          DPARAM(3) = -1.0
          SA7    B5+5
          EQ     B92           GO TO 92(B)
 B84      BX6    X2
          MX7    0
          SA6    B5+2          DPARAM(2) = 1.0
          SA7    B5+3
          SA6    B5+8          DPARAM(5) = 1.0
          SA7    B5+9
          BX6    -X6
 B92      SA7    B5+1          DPARAM(1) = -1.0
          SA6    B5
 B96      SA5    BIGINV        (X5,) = BIGINV
          FX1    X4*X5         DD1 = DD1 / (SQRBIG*SQRBIG)
          DX7    X3*X5
          FX6    X3*X5
          FX1    X1+X7
          DX7    X1+X6
          FX6    X1+X6         (X6,X7) = (X3,X4) * (X5,0)
          SA7    B1+1          DD1 = (X6,X7)
          SA6    B1
          SA2    B3            (X2,X3) = DB1
          SA3    B3+1
          SA4    SQRBIG        (X4,) = SQRBIG
          FX1    X3*X4         (X6,X7) = DB1*SQRBIG
          DX7    X2*X4
          FX6    X2*X4
          FX1    X1+X7
          DX7    X1+X6
          FX6    X1+X6         (X6,X7) = (X2,X3) * (X4,0)
          SA7    B3+1          DB1 = (X6,X7)
          SA6    B3
          SA2    B5+2          (X2,X3) = DPARAM(2)
          SA3    B5+3
          FX1    X3*X4         (X6,X7) = DPARAM(2)*SQRBIG
          DX7    X2*X4
          FX6    X2*X4
          FX1    X1+X7
          DX7    X1+X6
          FX6    X1+X6         (X6,X7) = (X2,X3) * (X4,0)
          SA7    B5+3          DPARAM(2) = (X6,X7)
          SA6    B5+2
          SA2    B5+6          (X2,X3) = DPARAM(4)
          SA3    B5+7
          FX1    X3*X4         (X6,X7) = DPARAM(4)*SQRBIG
          DX7    X2*X4
          FX6    X2*X4
          FX1    X1+X7
          DX7    X1+X6
          FX6    X1+X6         (X6,X7) = (X2,X3) * (X4,0)
          SA7    B5+7          DPARAM(4) = (X6,X7)
          SA6    B5+6
          EQ     THIRTY6       GO TO 36
 FOURTY8  SA3    B2            (X3,X4) = DD2
          SA4    B2+1
          SA5    BIGINV        (X5,) = BIGINV
          BX0    X3
          AX0    59
          BX0    X0-X3         (X0,) = DABS(DD2)
          FX0    X5-X0         IF ( DABS(DD2) .GT. BIGINV ) GO TO 60
          NX0    X0
          NG     X0,SIXTY
          ZR     X3,OUT        IF(DD2 .EQ. 0.0) GO TO OUT
          SA1    B5            (X1,) = DPARAM(1)
          SA2    UNIT
          ZR     X1,C84        IF (DPARAM(1)) 96,84,88(C)
          NG     X1,C96
          BX6    X2                                        IC88
          MX7    0
          SA6    B5+6          DPARAM(4) = 1.0
          SA7    B5+7
          BX6    -X6
          SA6    B5+4          DPARAM(3) = -1.0
          SA7    B5+5
          EQ     C92           GO TO 92(C)
 C84      BX6    X2
          MX7    0
          SA6    B5+2          DPARAM(2) = 1.0
          SA7    B5+3
          SA6    B5+8          DPARAM(5) = 1.0
          SA7    B5+9
          BX6    -X6
 C92      SA7    B5+1          DPARAM(1) = -1.0
          SA6    B5
 C96      SA5    BIG           (X5,) = BIG
          FX1    X4*X5         DD2 = DD2 * (SQRBIG*SQRBIG)
          DX7    X3*X5
          FX6    X3*X5
          FX1    X1+X7
          DX7    X1+X6
          FX6    X1+X6         (X6,X7) = (X3,X4) * (X5,0)
          SA7    B2+1          DD2 = (X6,X7)
          SA6    B2
          SA2    B5+4          (X2,X3) = DPARAM(3)
          SA3    B5+5
          SA4    SQRBIGI       (X4,) = SQRBIGI
          FX1    X3*X4         (X6,X7) = DPARAM(3)/SQRBIG
          DX7    X2*X4
          FX6    X2*X4
          FX1    X1+X7
          DX7    X1+X6
          FX6    X1+X6         (X6,X7) = (X2,X3) * (X4,0)
          SA7    B5+5          DPARAM(3) = (X6,X7)
          SA6    B5+4
          SA2    B5+8          (X2,X3) = DPARAM(5)
          SA3    B5+9
          FX1    X3*X4         (X6,X7) = DPARAM(5)/SQRBIG
          DX7    X2*X4
          FX6    X2*X4
          FX1    X1+X7
          DX7    X1+X6
          FX6    X1+X6         (X6,X7) = (X2,X3) * (X4,0)
          SA7    B5+9          DPARAM(5) = (X6,X7)
          SA6    B5+8
          EQ     FOURTY8       GO TO 48
 SIXTY    SA3    B2            (X3,X4) = DD2
          SA4    B2+1
          SA5    BIG           (X5,) = BIG
          BX0    X3
          AX0    59
          BX0    X0-X3         (X0,) = DABS(DD2)
          FX0    X0-X5         IF ( DABS(DD2) .LT. BIG ) RETURN
          NX0    X0
          NG     X0,OUT        GO TO OUT
          SA1    B5            (X1,) = DPARAM(1)
          SA2    UNIT
          ZR     X1,D84        IF (DPARAM(1)) 96,84,88(D)
          NG     X1,D96
          BX6    X2                                        ID88
          MX7    0
          SA6    B5+6          DPARAM(4) = 1.0
          SA7    B5+7
          BX6    -X6
          SA6    B5+4          DPARAM(3) = -1.0
          SA7    B5+5
          EQ     D92           GO TO 92(D)
 D84      BX6    X2
          MX7    0
          SA6    B5+2          DPARAM(2) = 1.0
          SA7    B5+3
          SA6    B5+8          DPARAM(5) = 1.0
          SA7    B5+9
          BX6    -X6
 D92      SA7    B5+1          DPARAM(1) = -1.0
          SA6    B5
 D96      SA5    BIGINV        (X5,) = BIGINV
          FX1    X4*X5         DD2 = DD2 / (SQRBIG*SQRBIG)
          DX7    X3*X5
          FX6    X3*X5
          FX1    X1+X7
          DX7    X1+X6
          FX6    X1+X6         (X6,X7) = (X3,X4) * (X5,0)
          SA7    B2+1          DD2 = (X6,X7)
          SA6    B2
          SA2    B5+4          (X2,X3) = DPARAM(3)
          SA3    B5+5
          SA4    SQRBIG        (X4,) = SQRBIG
          FX1    X3*X4         (X6,X7) = DPARAM(3)*SQRBIG
          DX7    X2*X4
          FX6    X2*X4
          FX1    X1+X7
          DX7    X1+X6
          FX6    X1+X6         (X6,X7) = (X2,X3) * (X4,0)
          SA7    B5+5          DPARAM(3) = (X6,X7)
          SA6    B5+4
          SA2    B5+8          (X2,X3) = DPARAM(5)
          SA3    B5+9
          FX1    X3*X4         (X6,X7) = DPARAM(5)*SQRBIG
          DX7    X2*X4
          FX6    X2*X4
          FX1    X1+X7
          DX7    X1+X6
          FX6    X1+X6         (X6,X7) = (X2,X3) * (X4,0)
          SA7    B5+9          DPARAM(5) = (X6,X7)
          SA6    B5+8
          EQ     SIXTY         GO TO 60
 OUT      OUTFTN DROTMG        RETURN
*
 P1       BSS    2
 P2       BSS    2
 TEMP     BSS    2
 BIG      DATA   17504000000000000000B
 BIGINV   DATA   16704000000000000000B
 RTWO     DATA   17214000000000000000B
 SQRBIG   DATA   17344000000000000000B
 SQRBIGI  DATA   17044000000000000000B
 TOL      DATA   0.0
 UNIT     DATA   17204000000000000000B
*
          END
*DECK,SROTM
          IDENT  SROTM
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL SROTM(N,SX,INCX,SY,INCY,SPARAM)
*
*         APPLY THE TWO-MULTIPLY,TWO-ADD,GIVENS TRANSFORMATION
*
*         TO 2XN MATRIX  (SXI1  ... SXIN )
*                        (SYI1  ... SYIN )
*
*         SXII  = SX(1 + (I-1)*INCX)  IF INCX .GE. 0
*               = SX(1 + (I-N)*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR SYII
*
*         CONTENTS OF SPARAM( ) MUST BE PREVIOUSLY DEFINED BY
*         SROTMG
*
*         SPARAM(1) = FLAG , INDICATES THE FORM OF THE MATRIX H
*         SPARAM(2) = H11
*         SPARAM(3) = H21
*         SPARAM(4) = H12
*         SPARAM(5) = H22
*
*         THE FLAG VALUES AND THE CORRESPONDING FORMS OF THE MATRIX H
*         -2. (1 0)   -1. (H11 H12)   0. ( 1 H12)   1. (H11 1 )
*             (0 1)       (H21 H22)      (H21 1 )      (-1 H22)
*
*
*         SX( ),SY( )               SINGLE PRECISION
*         N,INCX,INCY               INTEGER TYPE
*         SPARAM( )                 SINGLE PRECISION
*
*         ROUNDED ARITHMETIC INSTRUCTIONS ARE USED
*
*         WRITTEN BY  RICHARD J. HANSON
*                     SANDIA LABORATORIES
*                     ALBUQUERQUE, NEW MEXICO
***       1 JUNE 77
*
          ENTRY  SROTM
          VFD    42/5HSROTM,18/6
*
 SROTM    DATA   0
          INFTN  SROTM,6     PROPER LINKAGE (RUN,FTN) MACRO
*
          SA1    B1          (X1)=N
          SB7    1           (B7)=1
*
          SB1    X1          (B1)=N
          SB1    B1-B7       (B1)=N-1
*
          MI     B1,OUT      IF (N .LE. 0), QUIT
*
*
          SA1    B2          (X1)=SX(1)
          SA2    B3          (X2)=INCX
*
          SA3    B4          (X3)=SY(1)
          SA4    B5          (X4)=INCY
          SA5    B6          (X5)=SPARAM(1), (A5)=LOC(SPARAM(1))
*
          ZR     B1,INCYNN   IF (N .EQ. 1) NO NEED TO TEST FOR NEG. INC.
          SX0    -B1         (X0)=-(N-1)
*
*
          SB2    X2          (B2)=INCX
          SB3    X4          (B3)=INCY
*
          GE     B2,INCXNN   IF (INCX .GE. 0) NO ADDRESS FIXUP NEEDED
          DX2    X0*X2       COMPUTE -(N-1)*INCX
          SB7    A1          (B7)=LOC(SX(1))
          SA1    B7+X2       (X1)=SX(1+(1-N)*INCX),(A1)=LOC(X(1))
*
 INCXNN   GE     B3,INCYNN   IF (INCY .GE. 0) NO ADDRESS FIXUP NEEDED
          DX4    X0*X4       COMPUTE -(N-1)*INCY
          SB7    A3          (B7)=LOC(SY(1))
          SA3    B7+X4       (X3)=SY(1+(1-N)*INCY),(A3)=LOC(Y(1))
*
 INCYNN   SB7    1           (B7)=1
          SB6    B7          (B6)=I=1
          ZR     X5,SP1E0    IF (SPARAM(1) .EQ. 0.0)
          PL     X5,SP1E1    IF (SPARAM(1) .EQ. 1.0)
          SA4    STWO        (X4)=2.0
*
          RX4    X4+X5       (X4)=SPARAM(1)+2.0
          NX4    X4          (X4)=NORM.(X4)
          ZR     X4,OUT      IF (SPARAM(1) .EQ. -2.0), QUIT
*
*    HERE SPARAM(1)=-1.0.  PERFORM (RARELY USED) RESCALING LOOP
          SA2    A5+1        (X2)=SPARAM(2)=H11
          SA4    A5+3        (X4)=SPARAM(4)=H12
          BX0    X2          (X0)=H11
          SA2    A5+2        (X2)=SPARAM(3)=H21
          SA5    A5+4        (X5)=SPARAM(5)=H22
*
*    APPLY  (H11   H12)  TO (SX(1) ... SX(N))
*           (         )     (               )
*           (H21   H22)     (SY(1) ... SY(N))
          GT     B6,B1,CLR   IF (I .GT. N-1) CLEAN-UP LOGIC
 LOOP     RX6    X0*X1       (X6)=H11*SX(I)
          RX7    X1*X2       (X7)=H21*SX(I)
          RX1    X3*X4       (X1)=H12*SY(I)
          RX3    X3*X5       (X3)=H22*SY(I)
          RX6    X1+X6       (X6)=H11*SX(I)+H12*SY(I)
*
          SA1    A1+B2       (X1)=SX(I+1). NEXT ITER.
          NX6    X6          (X6)=NORM.(X6)
          RX7    X3+X7       (X7)=H21*SX(I)+H22*SY(I)
          SA3    A3+B3       (X3)=SY(I+1). NEXT ITER.
          SA6    A1-B2       SX(I)=(X6)
          NX7    X7          (X7)=NORM.(X7)
          SB6    B6+B7       (B6)=I=I+1. INCREMENT I.
          SA7    A3-B3       SY(I-1)=(X7)
*
          LE     B6,B1,LOOP  IF (I .LE. N-1) CONTINUE LOOP
 CLR      RX6    X0*X1       (X6)=H11*SX(N)
          RX7    X1*X2       (X7)=H21*SX(N)
          RX1    X3*X4       (X1)=H12*SY(N)
          RX3    X3*X5       (X3)=H22*SY(N)
          RX6    X1+X6       (X6)=H11*SX(N)+H12*SY(N)
          RX7    X3+X7       (X7)=H21*SX(N)+H22*SY(N)
          NX6    X6          (X6)=NORM.(X6)
          NX7    X7          (X7)=NORM.(X7)
          SA6    A1          SX(N)=(X6)
          SA7    A3          SY(N)=(X7)
          JP     OUT         QUIT
*
*    APPLY  ( 1    H12)  TO (SX(1) ... SX(N))
*           (         )     (               )
*           (H21    1 )     (SY(1) ... SY(N))
 SP1E0    SA2    A5+2        (X2)=SPARAM(3)=H21
          SA5    A5+3        (X5)=SPARAM(4)=H12
          BX0    X2          (X0)=H21
          SB1    B1-B7       (B1)=N-2
          GT     B6,B1,FIXN0 IF (I .GT. N-2) CLEAN-UP LOGIC
*
 LOOP0    SA2    A1+B2       (X2)=SX(I+1)
          SA4    A3+B3       (X4)=SY(I+1)
          RX7    X0*X1       (X7)=H21*SX(I)
          RX6    X3*X5       (X6)=H12*SY(I)
          SB6    B6+B7       (B6)=I=I+1. INCREMENT I.
          NO     3           DEAD
*
          RX7    X3+X7       (X7)=SY(I-1)+H21*SX(I-1)
          RX3    X0*X2       (X3)=H21*SX(I)
          RX6    X1+X6       (X6)=SX(I-1)+H12*SY(I-1)
          RX1    X4*X5       (X1)=H12*SY(I)
          NO     0           DEAD
          NX7    X7          (X7)=NORM.(X7)
          RX4    X4+X3       (X2)=SY(I)+H21*SX(I)
          NX6    X6          (X6)=NORM.(X6)
*
          SA3    A4+B3       (X3)=SY(I+1) NEXT ITERATION
          SA7    A4-B3       SY(I-1)=(X7)
          RX2    X1+X2       (X4)=SX(I)+H12*SY(I)
          SA6    A2-B2       SX(I-1)=(X6)
          NX7    X4          (X7)=NORM.(X4)
          SA1    A2+B2       (X1)=SX(I+1) NEXT ITERATION
          NX6    X2          (X6)=NORM.(X2)
          NO     0           DEAD
          SA7    A4          SY(I)=(X7)
          SB6    B6+B7       (B6)=I=I+1. INCREMENT I.
          NO     0           DEAD
          SA6    A2          SX(I-1)=(X6)
          LE     B6,B1,LOOP0 IF (I .LE. N-2) CONTINUE LOOP
 FIXN0    SB1    B1+B7       (B1)=N-1
          SB1    B1+B7       (B1)=N
*    HERE ONE VECTOR IS PRE-FETCHED. AT MOST TWO COMPS. REMAIN
 CL0      RX7    X0*X1       (X7)=H21*SX(I)
          RX6    X3*X5       (X6)=H12*SY(I)
          SB6    B6+B7       (B6)=I=I+1. INCREMENT I.
          RX7    X3+X7       (X7)=SY(I-1)+H21*SX(I-1)
          RX6    X1+X6       (X6)=SX(I-1)+H12*SY(I-1)
          NX7    X7          (X7)=NORM.(X7)
          NX6    X6          (X6)=NORM.(X6)
          SA7    A3          SY(I-1)=(X7)
          SA6    A1          SX(I-1)=(X6)
          GT     B6,B1,OUT   IF (I .GT. N) QUIT
          SA1    A1+B2       (X1)=SX(I)
          SA3    A3+B3       (X3)=SY(I)
          JP     CL0
*
*    APPLY  (H11    1 )  TO (SX(1) ... SX(N))
*           (         )     (               )
*           (-1    H22)     (SY(1) ... SY(N))
 SP1E1    SA2    A5+1        (X2)=SPARAM(2)=H11
          SA5    A5+4        (X5)=SPARAM(5)=H22
          BX0    X2          (X0)=H11
          SB1    B1-B7       (B1)=N-2
          GT     B6,B1,FIXN1 IF (I .GT. N-2) CLEAN-UP LOGIC
*
 LOOP1    SA2    A1+B2       (X2)=SX(I+1)
          SA4    A3+B3       (X4)=SY(I+1)
          RX7    X3*X5       (X7)=H22*SY(I)
          RX6    X0*X1       (X6)=H11*SX(I)
          SB6    B6+B7       (B6)=I=I+1. INCREMENT I.
          NO     3           DEAD
*
          RX7    X7-X1       (X7)=-SX(I-1)+H22*SY(I-1)
          RX1    X0*X2       (X1)=H11*SX(I)
          RX6    X3+X6       (X6)=SY(I-1)+H11*SX(I-1)
          RX3    X4*X5       (X3)=H22*SY(I)
          NO     0           DEAD
          NX7    X7          (X7)=NORM.(X7)
          RX4    X1+X4       (X4)=SY(I)+H11*SX(I)
          NX6    X6          (X6)=NORM.(X6)
*
          SA7    A4-B3       SY(I-1)=(X7)
          RX2    X3-X2       (X2)=-SX(I)+H22*SY(I)
          SA3    A4+B3       (X3)=SY(I+1) NEXT ITERATION
          SA6    A2-B2       SX(I-1)=(X6)
          NX6    X4          (X6)=NORM.(X4)
          SA1    A2+B2       (X1)=SX(I+1) NEXT ITERATION
          NO     0           DEAD
          NX7    X2          (X7)=NORM.(X2)
          SA6    A2          SX(I)=(X6)
          SB6    B6+B7       (B6)=I=I+1. INCREMENT I.
          NO     0           DEAD
          SA7    A4          SY(I-1)=(X7)
          LE     B6,B1,LOOP1 IF (I .LE. N-2) CONTINUE LOOP
 FIXN1    SB1    B1+B7       (B1)=N-1
          SB1    B1+B7       (B1)=N
*    HERE ONE VECTOR IS PRE-FETCHED. AT MOST TWO COMPS. REMAIN
 CL1      RX7    X0*X1       (X7)=H11*SX(I)
          RX6    X3*X5       (X6)=H22*SY(I)
          SB6    B6+B7       (B6)=I=I+1. INCREMENT I.
          RX7    X3+X7       (X7)=SY(I-1)+H11*SX(I-1)
          RX6    X6-X1       (X6)=-SX(I-1)+H22*SY(I-1)
          NX7    X7          (X7)=NORM.(X7)
          NX6    X6          (X6)=NORM.(X6)
          SA7    A1          SX(I-1)=(X7)
          SA6    A3          SY(I-1)=(X6)
          GT     B6,B1,OUT   IF (I .GT. N), QUIT
          SA1    A1+B2       (X1)=SX(I)
          SA3    A3+B3       (X3)=SY(I)
          JP     CL1
*
 OUT      OUTFTN SROTM
 STWO     DATA   2.0
*         END    SROTM
          END
*DECK,DROTM
          IDENT  DROTM
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL DROTM(N,DX,INCX,DY,INCY,DPARAM)
*
*         APPLY THE TWO-MULTIPLY,TWO-ADD,GIVENS TRANSFORMATION
*
*         TO 2XN MATRIX  (DXI1  ... DXIN )
*                        (DYI1  ... DYIN )
*
*         DXII  = DX(1 + (I-1)*INCX)  IF INCX .GE. 0
*               = DX(1 + (I-N)*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR DYII
*
*         CONTENTS OF DPARAM( ) MUST BE PREVIOUSLY DEFINED BY
*         DROTMG
*
*         DPARAM(1) = FLAG , INDICATES THE FORM OF THE MATRIX H
*         DPARAM(2) = H11
*         DPARAM(3) = H21
*         DPARAM(4) = H12
*         DPARAM(5) = H22
*
*         THE FLAG VALUES AND THE CORRESPONDING FORMS OF THE MATRIX H
*         -2. (1 0)   -1. (H11 H12)   0. ( 1 H12)   1. (H11 1 )
*             (0 1)       (H21 H22)      (H21 1 )      (-1 H22)
*
*
*         DX( ),DY( )               DOUBLE PRECISION
*         N,INCX,INCY               INTEGER TYPE
*         DPARAM( )                 DOUBLE PRECISION
*
*         WRITTEN BY  DAVID R. KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  DROTM
          VFD    42/5HDROTM,18/6
*
 DROTM    DATA   0             ENTRY/EXIT
          INFTN  DROTM,6
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          SB1    X1+B7         (B1) = N-1
*
          SA3    B3            (X3) = INCX
          NG     B1,OUT        IF N .LE. 0 , GO TO OUT
          SA5    B5            (X5) = INCY
          SX1    -B1           (X1) = -(N-1)
          LX3    1             INCX = 2*INCY
          IX5    X5+X5         INCY = 2*INCY
          SB3    X3            (B3) = INCX
          SB5    X5            (B5) = INCY
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3         LOC(DXI1 ) = LOC(DX)-(N-1)*INCX
          SB2    X3+B2         (B2) = LOC(DXI1 )
*
 ONE      GT     B5,TWO        IF INCY .GT. 0 ,GO TO TWO
          DX5    X1*X5         LOC(DYI1 ) = LOC(DY)-(N-1)*INCY
          SB4    X5+B4         (B4) = LOC(DYI1 )
*
 TWO      SA3    B6            (X3) = DPARAM(1)   IFLAG
          SA2    RTWO          (X2) = 2.0
          RX2    X3+X2
          NX2    X2
          ZR     X2,OUT        IF FLAG .EQ. -2.0 , GO TO OUT
*
          SA1    A3+2          (X1,X2) = DPARAM(2)
          SA2    A3+3
          BX6    X1            (X6,X7) = (X1,X2)
          BX7    X2
          SA6    H11           H11 = (X6,X7)
          SA7    H11+1
*
          SA1    A2+1          (X1,X2) = DPARAM(3)
          SA2    A2+2
          BX6    X1            (X6,X7) = (X1,X2)
          BX7    X2
          SA6    H21           H21 = (X6,X7)
          SA7    H21+1
*
          SA1    A2+1          (X1,X2) = DPARAM(4)
          SA2    A2+2
          BX6    X1            (X6,X7) = (X1,X2)
          BX7    X2
          SA6    H12           H12 = (X6,X7)
          SA7    H12+1
*
          SA1    A2+1          (X1,X2) = DPARAM(5)
          SA2    A2+2
          BX6    X1            (X6,X7) = (X1,X2)
          BX7    X2
          SA6    H22           H22 = (X6,X7)
          SA7    H22+1
*
          SB1    B1-B7         (B1) = N
          ZR     X3,LOOP1      IF FLAG .EQ. 0.0, GO TO LOOP1
          NG     X3,LOOP3      IF FLAG .EQ.-1.0, GO TO LOOP3
          EQ     LOOP2         IF FLAG .EQ. 1.0, GO TO LOOP2
*
*
 LOOP1    SA1    H12           (X1,X2) = H12      ITEN
          SA2    H12+1
          SA3    B4            (X3,X4) = DYII
          SA4    B4-B7
*
          FX5    X2*X3         (X0,X3) = H12*DYII
          FX0    X1*X4
          FX5    X0+X5
          FX4    X1*X3
          DX0    X1*X3
          FX5    X0+X5
          FX0    X4+X5
          DX3    X4+X5
*
          SA1    B2            (X1,X2) = DXII
          SA2    B2-B7
          BX6    X1            (X6,X7) = DXII
          BX7    X2
*
          FX4    X6+X0         (X6,X7) = (X6,X7)+(X0,X3)
          DX5    X6+X0
          FX0    X7+X3
          NX4    X4
          FX3    X0+X5
          FX0    X3+X4
          NX5    X0
          DX3    X3+X4
          NX4    X3
          FX6    X4+X5
          DX7    X4+X5
*
          SA6    A1            DXII  = (X6,X7)
          SA7    A2
*
          SA3    H21           (X3,X4) = H21
          SA4    H21+1
*
          FX5    X2*X3         (X6,X7) = H21*(X1,X2)
          FX0    X1*X4
          FX5    X0+X5
          FX4    X1*X3
          DX0    X1*X3
          FX5    X0+X5
          FX6    X4+X5
          DX7    X4+X5
*
          SA1    B4            (X1,X2) = DYII
          SA2    B4-B7
*
          FX4    X6+X1         (X6,X7) = (X6,X7)+(X1,X2)
          DX5    X6+X1
          FX1    X7+X2
          NX4    X4
          FX2    X1+X5
          FX1    X2+X4
          NX5    X1
          DX2    X2+X4
          NX4    X2
          FX6    X4+X5
          DX7    X4+X5
*
          SA6    A1            DYII  = (X6,X7)
          SA7    A2
*
          SB2    B2+B3         (B2) = LOC(DXII+1 )
          SB4    B4+B5         (B4) = LOC(DYII+1 )
*
          SB1    B1+B7         COUNT TERM
          NZ     B1,LOOP1      IF I .NE. N , LOOP1
*
          EQ     OUT           GO TO OUT
*
*
*
 LOOP2    SA1    B2            (X1,X2) = DXII      ITHIRTY
          SA2    B2-B7
          SA3    H11           (X3,X4) = H11
          SA4    H11+1
*
          FX5    X2*X3         (X0,X3) = H11*DXII
          FX0    X1*X4
          FX5    X0+X5
          FX4    X1*X3
          DX0    X1*X3
          FX5    X0+X5
          FX0    X4+X5
          DX3    X4+X5
*
          SA4    B4            (X4,X5) = DYII
          SA5    B4-B7
          BX6    X4
          BX7    X5            (X6,X7) = DYII
*
          FX4    X6+X0         (X6,X7) = (X6,X7)+(X0,X3)
          DX5    X6+X0
          FX0    X7+X3
          NX4    X4
          FX3    X0+X5
          FX0    X3+X4
          NX5    X0
          DX3    X3+X4
          NX4    X3
          FX6    X4+X5
          DX7    X4+X5
*
          SA6    A1            DXII  = (X6,X7)
          SA7    A2
*
          SA3    H22           (X3,X4) = H22
          SA4    H22+1
*
          BX6    X1            (X6,X7) = (X1,X2)    ISAVE OLD DX
          BX7    X2
*
          SA1    B4            (X1,X2) = DYII
          SA2    B4-B7
*
          FX5    X2*X3         (X1,X2) = H22*DYII
          FX0    X1*X4
          FX5    X0+X5
          FX4    X1*X3
          DX0    X1*X3
          FX5    X0+X5
          FX1    X4+X5
          DX2    X4+X5
*
          FX4    X1-X6         (X6,X7) = -(X6,X7)+(X1,X2)
          DX5    X1-X6
          FX1    X2-X7
          NX4    X4
          FX2    X1+X5
          FX1    X2+X4
          NX5    X1
          DX2    X2+X4
          NX4    X2
          FX6    X4+X5
          DX7    X4+X5
*
          SA6    A1            DYII  = (X6,X7)
          SA7    A2
*
          SB2    B2+B3         (B2) = LOC(DXII+1 )
          SB4    B4+B5         (B4) = LOC(DYII+1 )
*
          SB1    B1+B7         COUNT TERM
          NZ     B1,LOOP2      IF I .NE. N , LOOP2
*
          EQ     OUT           GO TO OUT
*
*
*
 LOOP3    SA1    B2            (X1,X2) = DXII     IFIFTY
          SA2    B2-B7
          SA3    H11           (X3,X4) = H11
          SA4    H11+1
*
          FX5    X2*X3         (X6,X7) = DXII *H11
          FX0    X1*X4
          FX5    X0+X5
          FX4    X1*X3
          DX0    X1*X3
          FX5    X0+X5
          FX6    X4+X5
          DX7    X4+X5
*
          SA1    B4            (X1,X2) = DYII
          SA2    B4-B7
          SA3    H12           (X3,X4) = H12
          SA4    H12+1
*
          FX5    X2*X3         (X1,X2) = DYII *H12
          FX0    X1*X4
          FX5    X0+X5
          FX4    X1*X3
          DX0    X1*X3
          FX5    X0+X5
          FX1    X4+X5
          DX2    X4+X5
*
          FX4    X6+X1         (X6,X7) = (X6,X7)+(X1,X2)
          DX5    X6+X1
          FX1    X7+X2
          NX4    X4
          FX2    X1+X5
          FX1    X2+X4
          NX5    X1
          DX2    X2+X4
          NX4    X2
          FX6    X4+X5
          DX7    X4+X5
*
          SA6    DW            DW = (X6,X7)
          SA7    DW+1
*
          SA1    B2            (X1,X2) = DXII
          SA2    B2-B7
*
          SA3    H21           (X3,X4) = H21
          SA4    H21+1
*
          FX5    X2*X3         (X6,X7) = DXII *H21
          FX0    X1*X4
          FX5    X0+X5
          FX4    X1*X3
          DX0    X1*X3
          FX5    X0+X5
          FX6    X4+X5
          DX7    X4+X5
*
          SA1    B4            (X1,X2) = DYII
          SA2    B4-B7
          SA3    H22           (X3,X4) = H22
          SA4    H22+1
*
          FX5    X2*X3         (X1,X2) = DYII *H22
          FX0    X1*X4
          FX5    X0+X5
          FX4    X1*X3
          DX0    X1*X3
          FX5    X0+X5
          FX1    X4+X5
          DX2    X4+X5
*
*
          FX4    X6+X1         (X6,X7) = (X6,X7)+(X1,X2)
          DX5    X6+X1
          FX1    X7+X2
          NX4    X4
          FX2    X1+X5
          FX1    X2+X4
          NX5    X1
          DX2    X2+X4
          NX4    X2
          FX6    X4+X5
          DX7    X4+X5
*
          SA6    A1            DYII  = (X6,X7)
          SA7    A2
*
          SA3    DW            (X3,X4) = DW
          SA4    DW+1
          BX6    X3            (X6,X7) = (X3,X4)
          BX7    X4
          SA6    B2            DXII  = (X6,X7)
          SA7    B2+1
*
          SB1    B1+B7         COUNT TERM
          SB2    B2+B3         (B2) = LOC(DXII+1 )
          SB4    B4+B5         (B4) = LOC(DYII+1 )
*
          NZ     B1,LOOP3      IF I .NE. N ,LOOP3
*
 OUT      OUTFTN DROTM         RETURN
*
 DW       BSS    2
 H11      BSS    2
 H21      BSS    2
 H12      BSS    2
 H22      BSS    2
*
 RTWO     DATA   2.0
*
          END
*DECK,SCOPY
          IDENT  SCOPY
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL SCOPY(N,SX,INCX,SY,INCY)
*
*         COPY VECTOR ELEMENT SXII  INTO SYII  FOR I=1 TO N
*
*         SXII  = SX(1 + (I-1)*INCX)  IF INCX .GE. 0
*               = SX(1 + (I-N)*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR SYII
*
*         SX( ),SY( )               SINGLE PRECISION
*         N,INCX,INCY               INTEGER TYPE
*
*         WRITTEN BY  DAVID R. KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  SCOPY
          VFD    42/5HSCOPY,18/5
*
 SCOPY    DATA   0             ENTRY/EXIT
          INFTN  SCOPY,5
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          SB1    X1+B7         (B1) = N-1
*
          SA3    B3            (X3) = INCX
          NG     B1,OUT        IF N .LE. O , GO TO OUT
          SA5    B5            (X5) = INCY
          SX1    -B1           (X1) = -(N-1)
          SB3    X3            (B3) = INCX
          SB5    X5            (B5) = INCY
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3         LOC(SXI1 ) = LOC(SX) - (N-1)*INCX
          SB2    X3+B2         (B2) = LOC(SXI1 )
*
 ONE      GT     B5,TWO        IF INCY .GT. 0 , GO TO TWO
          DX5    X1*X5         LOC(SYI1 ) = LOC(SY) - (N-1)*INCY
          SB4    X5+B4         (B4) = LOC(SYI1 )
*
*                              (I = 1)
 TWO      SA2    B2            (X2) = SXI1
          SA4    B4            (A4) = LOC(SYI1 )
          BX6    X2
          SA6    B4            SXI1  TO SYI1
          ZR     B1,OUT        IF I .EQ. N , GO TO OUT
*
*                              (I = I+1)
 LOOP     SA2    A2+B3         (X2) = SXII
          SA4    A4+B5         (A4) = LOC(SYII )
          BX6    X2
          SB1    B1+B7         COUNT TERM
          SA6    A4            SXII  TO SYII
          NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
 OUT      OUTFTN SCOPY         RETURN
          END
*DECK,DCOPY
          IDENT  DCOPY
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL DCOPY(N,DX,INCX,DY,INCY)
*
*         COPY VECTOR ELEMENT DXII  INTO DYII  FOR I=1 TO N
*
*         DXII  = DX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*               = DX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR DYII
*
*         DX( ),DY( )               DOUBLE PRECISION
*         N,INCX,INCY               INTEGER TYPE
*
*         WRITTEN BY  DAVID R.KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  DCOPY
          VFD    42/5HDCOPY,18/5
*
 DCOPY    DATA   0             ENTRY/EXIT
          INFTN  DCOPY,5
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          SB1    X1+B7         (B1) = N-1
*
          SA3    B3            (X3) = INCX
          NG     B1,OUT        IF N .LE. 0 , GO TO OUT
          SA5    B5            (X5) = INCY
          SX1    -B1           (X1) = -(N-1)
          LX3    1             INCX = 2*INCX
          IX5    X5+X5         INCY = 2*INCY
          SB3    X3            (B3) = INCX
          SB5    X5            (B5) = INCY
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3         LOC(DXI1 ) = LOC(DX) - (N-1)*INCX
          SB2    X3+B2         (B2) = LOC(DXI1 )
*
 ONE      GT     B5,TWO        IF INCY .GT. 0 , GO TO TWO
          DX5    X1*X5         LOC(DYI1 ) = LOC(DY) - (N-1)*INCY
          SB4    X5+B4         (B4) = LOC(DYI1 )
*
*                              (I = 1)
 TWO      SA2    B2            (X2) = DXI1
          SA4    B4            (A4) = LOC(DYI1 )
          BX6    X2
          SA5    B2-B7         (X4,X5) = DXI1
          SA6    B4
          BX7    X5
          SA7    B4-B7         DXI1  TO DYI1
          ZR     B1,OUT        IF I .EQ. N , GO TO OUT
*
*                              (I = I+1)
 LOOP     SA2    A2+B3         (X2) = DXII
          SA4    A4+B5         (A4) = LOC(DYII )
          BX6    X2
          SA5    A2-B7         (X4,X5) = DXII
          SA6    A4
          BX7    X5
          SB1    B1+B7         COUNT TERM
          SA7    A4-B7         DXII  TO DYII
          NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
 OUT      OUTFTN DCOPY         RETURN
          END
*DECK,CCOPY
          IDENT  CCOPY
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL CCOPY(N,CX,INCX,CY,INCY)
*
*         COPY VECTOR ELEMENT CXII  INTO CYII  FOR I=1 TO N
*
*         CXII  = CX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*               = CX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR CYII
*
*         CX( ),CY( )               COMPLEX TYPE
*         N,INCX,INCY               INTEGER TYPE
*
*         WRITTEN BY  DAVID R.KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  CCOPY
          VFD    42/5HCCOPY,18/5
*
 CCOPY    DATA   0             ENTRY/EXIT
          INFTN  CCOPY,5
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          SB1    X1+B7         (B1) = N-1
*
          SA3    B3            (X3) = INCX
          NG     B1,OUT        IF N .LE. 0 , GO TO OUT
          SA5    B5            (X5) = INCY
          SX1    -B1           (X1) = -(N-1)
          LX3    1             INCX = 2*INCX
          IX5    X5+X5         INCY = 2*INCY
          SB3    X3            (B3) = INCX
          SB5    X5            (B5) = INCY
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3         LOC(CXI1 ) = LOC(CX) - (N-1)*INCX
          SB2    X3+B2         (B2) = LOC(CXI1 )
*
 ONE      GT     B5,TWO        IF INCY .GT. 0 , GO TO TWO
          DX5    X1*X5         LOC(CYI1 ) = LOC(CY) - (N-1)*INCY
          SB4    X5+B4         (B4) = LOC(CYI1 )
*
*                              (I = 1)
 TWO      SA2    B2            (X2) = CXI1
          SA4    B4            (A4) = LOC(CYI1 )
          BX6    X2
          SA5    B2-B7         (X4,X5) = CXI1
          SA6    B4
          BX7    X5
          SA7    B4-B7         CXI1  TO CYII
          ZR     B1,OUT        IF I .EQ. N , GO TO OUT
*
*                              (I = I+1)
 LOOP     SA2    A2+B3         (X2) = CXII
          SA4    A4+B5         (A4) = LOC(CYII )
          BX6    X2
          SA5    A2-B7         (X4,X5) = CXII
          SA6    A4
          BX7    X5
          SB1    B1+B7         COUNT TERM
          SA7    A4-B7         CXII  TO CYII
          NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
 OUT      OUTFTN CCOPY         RETURN
          END
*DECK,SSWAP
          IDENT  SSWAP
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL SSWAP(N,SX,INCX,SY,INCY)
*
*         INTERCHANGE VECTOR ELEMENTS SXII  AND SYII  FOR I=1 TO N
*
*         SXII  = SX(1 + (I-1)*INCX)  IF INCX .GE. 0
*               = SX(1 + (I-N)*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR SYII
*
*         SX( ),SY( )               SINGLE PRECISION
*         N,INCX,INCY               INTEGER TYPE
*
*         WRITTEN BY  DAVID R. KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  SSWAP
          VFD    42/5HSSWAP,18/5
*
 SSWAP    DATA   0             ENTRY/EXIT
          INFTN  SSWAP,5
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          SB1    X1+B7         (B1) = N-1
*
          SA3    B3            (X3) = INCX
          NG     B1,OUT        IF N .LE. 0 , GO TO OUT
          SA5    B5            (X5) = INCY
          SX1    -B1           (X1) = -(N-1)
          SB3    X3            (B3) = INCX
          SB5    X5            (B5) = INCY
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3         LOC(XI1 ) = LOC(SX) - (N-1)*INCX
          SB2    X3+B2         (B2) = LOC(SXI1 )
*
 ONE      GT     B5,TWO        IF INCY .GT. 0 , GO TO TWO
          DX5    X1*X5         LOC(YI1 ) = LOC(SY) - (N-1)*INCY
          SB4    X5+B4         (B4) = LOC(SYI1 )
*
*                              (I = 1)
 TWO      SA2    B2            (X2) = SXI1
          SA4    B4            (X4) = SYI1
          BX6    X2            (X6) = (X2)
          BX7    X4            (X7) = (X4)
          SA6    B4            SXI1  TO SYI1
          SA7    B2            SYI1  TO SXI1
          ZR     B1,OUT        IF I .EQ. N , GO TO OUT
*
*                              (I = I+1)
 LOOP     SA2    A2+B3         (X2) = SXII
          SA4    A4+B5         (X4) = SYII
          BX6    X2            (X6) = (X2)
          BX7    X4            (X7) = (X4)
          SB1    B1+B7         COUNT TERM
          SA6    A4            SXII  TO SYII
          SA7    A2            SYII  TO SXII
          NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
 OUT      OUTFTN SSWAP         RETURN
          END
*DECK,DSWAP
          IDENT  DSWAP
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL DSWAP(N,DX,INCX,DY,INCY)
*
*         INTERCHANGE VECTOR ELEMENTS DXII  AND DYII  FOR I=1 TO N
*
*         DXII  = DX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*               = DX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR DYII
*
*         DX( ),DY( )               DOUBLE PRECISION
*         N,INCX,INCY               INTEGER TYPE
*
*         WRITTEN BY  DAVID R. KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  DSWAP
          VFD    42/5HDSWAP,18/5
*
 DSWAP    DATA   0             ENTRY/EXIT
          INFTN  DSWAP,5
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          SB1    X1+B7         (B1) = N-1
*
          SA3    B3            (X3) = INCX
          NG     B1,OUT        IF N .LE.0 , GO TO OUT
          SA5    B5            (X5) = INCY
          SX1    -B1           (X1) = -(N-1)
          LX3    1             INCX = 2*INCX
          IX5    X5+X5         INCY = 2*INCY
          SB3    X3            (B3) = INCX
          SB5    X5            (B5) = INCY
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3         LOC(XI1 ) = LOC(DX) - (N-1)*INCX
          SB2    X3+B2         (B2) = LOC(DXI1 )
*
 ONE      GT     B5,TWO        IF INCY .GT. 0 , GO TO TWO
          DX5    X1*X5         LOC(YI1 ) = LOC(DY) - (N-1)*INCY
          SB4    X5+B4         (B4) = LOC(DYI1 )
*
*                              (I = 1)
 TWO      SA2    B2
          SA4    B4
          BX6    X2
          BX7    X4
          SA6    B4
          SA7    B2
*
          SA3    A2-B7         (X2,X3) = DXI1
          SA5    A4-B7         (X4,X5) = DYI1
          BX6    X3
          BX7    X5
          SA6    A5            DXI1  = DYI1
          SA7    A3            DYI1  = DXI1
          ZR     B1,OUT        IF I .EQ. N , GO TO OUT
*
*                              (I = I+1)
 LOOP     SA2    A2+B3
          SA4    A4+B5
          BX6    X2
          BX7    X4
          SA6    A4
          SA7    A2
*
          SA3    A2-B7         (X2,X3) = DXII
          SA5    A4-B7         (X4,X5) = DYII
          BX6    X3
          BX7    X5
          SB1    B1+B7         COUNT TERM
          SA6    A5            DXII  = DYII
          SA7    A3            DYII  = DXII
          NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
 OUT      OUTFTN DSWAP         RETURN
          END
*DECK,CSWAP
          IDENT  CSWAP
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL CSWAP(N,CX,INCX,CY,INCY)
*
*         INTERCHANGE VECTOR ELEMENTS CXII  AND CYII  FOR I=1 TO N
*
*         CXII  = CX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*               = CX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         SIMILAR DEFINITIONS FOR CYII
*
*         CX( ),CY( )               COMPLEX TYPE
*         N,INCX,INCY               INTEGER TYPE
*
*         WRITTEN BY  DAVID R. KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  CSWAP
          VFD    42/5HCSWAP,18/5
*
 CSWAP    DATA   0             ENTRY/EXIT
          INFTN  CSWAP,5
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          SB1    X1+B7         (B1) = N-1
*
          SA3    B3            (X3) = INCX
          NG     B1,OUT        IF N .LE. N , GO TO OUT
          SA5    B5            (X5) = INCY
          SX1    -B1           (X1) = -(N-1)
          LX3    1             INCX = 2*INCX
          IX5    X5+X5         INCY = 2*INCY
          SB3    X3            (B3) = INCX
          SB5    X5            (B5) = INCY
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3         LOC(XI1 ) = LOC(CX) - (N-1)*INCX
          SB2    X3+B2         (B2) = LOC(CXI1 )
*
 ONE      GT     B5,TWO        IF INCY .GT. 0 , GO TO TWO
          DX5    X1*X5         LOC(YI1 ) = LOC(CY) - (N-1)*INCY
          SB4    X5+B4         (B4) = LOC(CYI1 )
*
*                              (I = 1)
 TWO      SA2    B2
          SA4    B4
          BX6    X2
          BX7    X4
          SA6    B4
          SA7    B2
*
          SA3    A2-B7         (X2,X3) = CXI1
          SA5    A4-B7         (X4,X5) = CYI1
          BX6    X3
          BX7    X5
          SA6    A5            CXI1  = CYI1
          SA7    A3            CYI1  = CXI1
          ZR     B1,OUT        IF I .EQ. N , GO TO OUT
*
*                              (I = I+1)
 LOOP     SA2    A2+B3
          SA4    A4+B5
          BX6    X2
          BX7    X4
          SA6    A4
          SA7    A2
*
          SA3    A2-B7         (X2,X3) = CXII
          SA5    A4-B7         (X4,X5) = CYII
          BX6    X3
          BX7    X5
          SB1    B1+B7         COUNT TERM
          SA6    A5            CXII  = CYII
          SA7    A3            CYII  = CXII
          NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
 OUT      OUTFTN CSWAP         RETURN
          END
*DECK,SNRM2
          IDENT  SNRM2
*
***       REAL FUNCTION  SNRM2(N,SX,INCX)
*
*         COMPUTES 2-VECTOR NORM (EUCLIDEAN NORM)
*
*         COMPUTED AS THE SQUARE ROOT OF THE SUM FROM I=1 TO N OF SXII *
*
*         SXII  = SX(1 + (I-1)*INCX)  IF INCX .GE. 0
*               = SX(1 + (I-N)*INCX)  IF INCX .LT. 0
*
*         SX( )                     SINGLE PRECISION
*         N,INCX                    INTEGER TYPE
*         SUM ACCUMULATED IN        SINGLE PRECISION
*         RESULT  SNRM2   IN        SINGLE PRECISION
*
*         ROUNDED ARITHMETIC INSTRUCTIONS ARE USED
*
*         WRITTEN BY  DAVID R. KINCAID AND ELIZABETH WILLIAMS
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  SNRM2
          VFD    42/5HSNRM2,18/3
*
 SNRM2    DATA   0             ENTRY/EXIT
          INFTN  SNRM2,3
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          MX6    0             (X6) = 0
          SB1    X1+B7         (B1) = N-1
*
          SA3    B3            (X3) = INCX
          NG     B1,OUT        IF N .LE. 0 , GO TO OUT
          SX1    -B1           (X1) = -(N-1)
          SB3    X3            (B3) = INCX
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3         LOC(SXI1 ) = LOC(SX) - (N-1)*INCX
          SB2    X3+B2         (B2) = LOC(SXI1 )
*
*                              (I = 1)
 ONE      SA2    B2            (X2) = SXI1
          RX1    X2*X2         (X2) = SXI1 *SXI1
*
          ZR     B1,EXIT       IF I .EQ. N , GO TO EXIT
*
*                              (I = I+1)
 LOOP     SA2    A2+B3         (X2) = SXII
          RX0    X1+X6         (X6) = (X6) + (X1)
          SB1    B1+B7         I = I+1
          NX6    X0
          RX1    X2*X2         (X1) = SXII *SXII
*
          NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
*                              (I = N)
 EXIT     RX0    X1+X6         (X6) = (X6) + (X1)
          NX6    X0
          SB1    RES           (B1) = LOC(RES)
          SA6    B1            RES  = (X6)
          CALL   SQRT,(B1)      (X6) =   SQRT(RES)
*
 OUT      OUTFTN SNRM2         RETURN
*
 RES      BSS    1
          END
*DECK,DNRM2
          IDENT  DNRM2
*
***       REAL FUNCTION  DNRM2(N,DX,INCX)
*
*         COMPUTES 2-VECTOR NORM (EUCLIDEAN NORM)
*
*         COMPUTED AS THE SQUARE ROOT OF THE SUM FROM I=1 TO N OF DXII *
*
*         DXII  = DX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*               = DX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         DX( )                     DOUBLE PRECISION
*         N,INCX                    INTEGER TYPE
*         SUM ACCUMULATED IN        DOUBLE PRECISION
*         RESULT  DNRM2   IN        DOUBLE PRECISION
*
*         WRITTEN BY  DAVID R. KINCAID AND ELIZABETH WILLIAMS
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  DNRM2
          VFD    42/5HDNRM2,18/3
*
 DNRM2    DATA   0             ENTRY/EXIT
          INFTN  DNRM2,3
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          MX6    0
          SB1    X1+B7         (B1) = N-1
          MX7    0             (X6,X7) = 0
*
          SA3    B3            (X3) = INCX
          NG     B1,OUT        IF N .LE. 0 , GO TO OUT
          SX1    -B1           (X1) = -(N-1)
          LX3    1             INCX = 2*INCX
          SB3    X3            (B3) = INCX
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3         LOC(DXI1 ) = LOC(DX) - (N-1)*INCX
          SB2    X3+B2         (B2) = LOC(DXI1 )
*
*
 ONE      SA1    B2            (X1,X2) = DXI1
          SA2    B2-B7
*
          FX0    X1*X2         (X0,X2) = DXI1 *DXI1
          FX5    X0+X0
          FX4    X1*X1
          DX0    X1*X1
          FX5    X0+X5
          FX0    X4+X5
          DX2    X4+X5
*
          ZR     B1,EXIT       IF I .EQ. N , GO TO EXIT
*
*                              (I = I+1)
 LOOP     SA1    A1+B3
*
          FX4    X6+X0         (X6,X7) = (X6,X7) + (X0,X2)
          DX5    X6+X0
          FX0    X7+X2
          NX4    X4
          FX2    X0+X5
          FX0    X2+X4
          NX5    X0
          DX2    X2+X4
          NX4    X2
          FX6    X4+X5
          DX7    X4+X5
*
          SA2    A1-B7         (X1,X2) = DXII
          SB1    B1+B7         I = I+1
*
          FX0    X1*X2         (X0,X2) = DXII *DXII
          FX5    X0+X0
          FX4    X1*X1
          DX0    X1*X1
          FX5    X0+X5
          FX0    X4+X5
          DX2    X4+X5
*
          NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
*                              (I = N)
 EXIT     FX4    X6+X0         (X6,X7) = (X6,X7) + (X0,X2)
          DX5    X6+X0
          FX0    X7+X2
          NX4    X4
          FX2    X0+X5
          FX0    X2+X4
          NX5    X0
          DX2    X2+X4
          NX4    X2
          FX6    X4+X5
          DX7    X4+X5
*
          SB1    RES           (B1) = RES
          SA6    B1            (RES) = (X6,X7)
          SA7    B1-B7
*
          CALL   DSQRT,(B1)     (X6,X7) =   SQRT(RES)
*
 OUT      OUTFTN DNRM2         RETURN
*
 RES      BSS    2
          END
*DECK,SCNRM2
          IDENT  SCNRM2
*
***       REAL FUNCTION  SCNRM2(N,CX,INCX)
*
*         COMPUTES 2-VECTOR NORM (EUCLIDEAN NORM)
*
*         COMPUTED AS THE SQUARE ROOT OF THE SUM
*         FROM I=1 TO N OF CONJ(CXII ) * CXII
*
*         CXII  = CX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*               = CX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         CX( )                     COMPLEX TYPE
*         N,INCX                    INTEGER TYPE
*         SUM ACCUMULATED IN        SINGLE PRECISION
*         RESULT  SCNRM2  IN        SINGLE PRECISION
*
*         ROUNDED ARITHMETIC INSTRUCTIONS ARE USED
*
*         WRITTEN BY  DAVID R. KINCAID AND ELIZABETH WILLIAMS
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  SCNRM2
          VFD    42/6HSCNRM2,18/3
*
 SCNRM2   DATA   0             ENTRY/EXIT
          INFTN  SCNRM2,3
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          MX6    0
          SB1    X1+B7         (B1) = N-1
*
          SA3    B3            (X3) = INCX
          NG     B1,OUT        IF N .LE. 0 , GO TO OUT
          SX1    -B1           (X1) = -(N-1)
          LX3    1             INCX = 2*INCX
          SB3    X3            (B3) = INCX
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          ZR     B3,OUT        IF INCX .EQ. 0 ,GO TO OUT
          DX3    X1*X3         LOC(CXI1 ) = LOC(CX) - (N-1)*INCX
          SB2    X3+B2         (B2) = LOC(CXI1 )
*
*                              (I = 1)
 ONE      SA1    B2            (X1) = REAL(CXI1 )
          SA2    B2-B7         (X2) = IMAG(CXI1 )
*
          RX0    X1*X1         (X0) = (REAL(CXI1 )**2
          RX5    X2*X2         (X5) = (IMAG(CXI1 )**2
          RX4    X0+X5         (X4) = (X0) + (X5)
          NX4    X4
*
          ZR     B1,EXIT       IF I .EQ. N , GO TO EXIT
*
*                              (I = I+1)
 LOOP     SA1    A1+B3         (X1) = REAL(CXII )
*
          RX5    X6+X4         (X6) = (X6) + (X4)
          SA2    A1-B7         (X2) = IMAG(CXII )
          NX6    X5
          RX0    X1*X1         (X0) = (REAL(CXII )**2
          RX5    X2*X2         (X5) = (IMAG(CXII )**2
          SB1    B1+B7         I = I+1
          RX4    X0+X5         (X4) = (X0) + (X5)
          NX4    X4
*
          NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
*                              (I = N)
 EXIT     RX5    X6+X4         (X6) = (X6) + (X4)
          NX6    X5
*
          SB1    RES           (B1) = RES
          SA6    B1            (RES) = (X6)
*
          CALL   SQRT,(B1)      (X6) =   SQRT(RES)
*
 OUT      OUTFTN SCNRM2        RETURN
*
 RES      BSS    1
          END
*DECK,SASUM
          IDENT  SASUM
*
***       REAL FUNCTION  SASUM(N,SX,INCX)
*
*         COMPUTES 1-VECTOR NORM
*
*         COMPUTED AS THE SUM FROM I=1 TO N OF THE ABSOLUTE VALUE OF SXI
*
*         SXII  = SX(1 + (I-1)*INCX)  IF INCX .GE. 0
*               = SX(1 + (I-N)*INCX)  IF INCX .LT. 0
*
*         SX( )                     SINGLE PRECISION
*         N,INCX                    INTEGER TYPE
*         SUM ACCUMULATED IN        SINGLE PRECISION
*         RESULT  SASUM   IN        SINGLE PRECISION
*
*         ROUNDED ARITHMETIC INSTRUCTIONS ARE USED
*
*         WRITTEN BY  DAVID R. KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  SASUM
          VFD    42/5HSASUM,18/3
*
 SASUM    DATA   0             ENTRY/EXIT
          INFTN  SASUM,3
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          MX6    0             (X6) = 0
          SB1    X1+B7         (B1) = N-1
*
          SA3    B3            (X3) = INCX
          NG     B1,OUT        IF N .LE. 0 , GO TO OUT
          SX1    -B1           (X1) = -(N-1)
          SB3    X3            (B3) = INCX
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3         LOC(RXI1 ) = LOC(RX) - (N-1)*INCX
          SB2    X3+B2         (B2) = LOC(RXI1 )
*
*                              (I=1)
 ONE      SA2    B2            (X2) = RXI1
          BX4    X2
          AX2    59
          BX5    X2-X4         (X5) = ABS(RXI1 )
*
          FX3    X6+X5         (X6) = (X6) + (X5)
          NX6    X3
*
          ZR     B1,OUT        IF I .EQ. N , GO TO OUT
*
*                              (I = I+1)
 LOOP     SA2    A2+B3         (X2) = RXII
          BX4    X2
          AX2    59
          BX5    X2-X4         (X5) = ABS(RXII )
*
          FX3    X6+X5         (X6) = (X6) + (X5)
          SB1    B1+B7         COUNT TERM
          NX6    X3
*
          NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
 OUT      OUTFTN SASUM         RETURN
          END
*DECK,DASUM
          IDENT  DASUM
*
***       REAL FUNCTION  DASUM(N,DX,INCX)
*
*         COMPUTES 1-VECTOR NORM
*
*         COMPUTED AS THE SUM FROM I=1 TO N OF THE ABSOLUTE VALUE OF DXI
*
*         DXII  = DX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*               = DX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         DX( )                     DOUBLE PRECISION
*         N,INCX                    INTEGER TYPE
*         SUM ACCUMULATED IN        DOUBLE PRECISION
*         RESULT  DASUM   IN        DOUBLE PRECISION
*
*         WRITTEN BY  DAVID R. KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  DASUM
          VFD    42/5HDASUM,18/3
*
 DASUM    DATA   0             ENTRY/EXIT
          INFTN  DASUM,3
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          MX6    0
          SB1    X1+B7         (B1) = N-1
          MX7    0             (X6,X7) = 0
*
          SA3    B3            (X3) = INCX
          NG     B1,OUT        IF N .LE. 0 , GO TO OUT
          LX3    1             INCX = 2*INCX
          SX1    -B1           (X1) = -(N-1)
          SB3    X3            (B3) = INCX
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3         LOC(DXI1 ) = LOC(DX) - (N-1)*INCX
          SB2    X3+B2         (B2) = LOC(DXI1 )
*                              (I=1)
 ONE      SA2    B2
          SA3    B2-B7         (X2,X3) = DXI1
          BX0    X2
          BX1    X3
          AX2    59
          AX3    59
          BX4    X2-X0
          BX5    X3-X1         (X4,X5) = DABS(DXI1 )
*
          FX0    X6+X4         (X6,X7) = (X6,X7) + (X4,X5)
          DX1    X6+X4
          FX4    X7+X5
          NX0    X0
          FX5    X4+X1
          FX4    X5+X0
          NX1    X4
          DX5    X5+X0
          NX0    X5
          FX6    X0+X1
          DX7    X0+X1
*
          ZR     B1,OUT        IF I .EQ. N , GO TO OUT
*
*                              (I = I+1)
 LOOP     SA2    A2+B3
          SA3    A2-B7         (X2,X3) = DXII
          BX0    X2
          BX1    X3
          AX2    59
          AX3    59
          BX4    X2-X0
          BX5    X3-X1         (X4,X5) = DABS(DXII )
*
          SB1    B1+B7         COUNT TERM
*
          FX0    X6+X4         (X6,X7) = (X6,X7) + (X4,X5)
          DX1    X6+X4
          FX4    X7+X5
          NX0    X0
          FX5    X4+X1
          FX4    X5+X0
          NX1    X4
          DX5    X5+X0
          NX0    X5
          FX6    X0+X1
          DX7    X0+X1
*
          NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
 OUT      OUTFTN DASUM         RETURN
          END
*DECK,SCASUM
          IDENT  SCASUM
*
***       REAL FUNCTION  SCASUM(N,CX,INCX)
*
*         COMPUTED AS THE SUM FROM I=1 TO N OF THE ABSOLUTE VALUE
*         OF REAL(CXII ) AND THE ABSOLUTE VALUE OF IMAG(CXII )
*
*         CXII  = CX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*               = CX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         CX( )                     COMPLEX TYPE
*         N,INCX                    INTEGER TYPE
*         SUM ACCUMULATED IN        SINGLE PRECISION
*         RESULT  SCASUM  IN        SINGLE PRECISION
*
*         ROUNDED ARITHMETIC INSTRUCTIONS ARE USED
*
*         WRITTEN BY  DAVID R. KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  SCASUM
          VFD    42/6HSCASUM,18/3
*
 SCASUM   DATA   0             ENTRY/EXIT
          INFTN  SCASUM,3
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          MX6    0
          SB1    X1+B7         (B1) = N-1
*
          SA3    B3            (X3) = INCX
          NG     B1,OUT        IF N .LE. 0 , GO TO OUT
          LX3    1             INCX = 2*INCX
          SX1    -B1           (X1) = -(N-1)
          SB3    X3            (B3) = INCX
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3         LOC(CXI1 ) = LOC(CX) - (N-1)*INCX
          SB2    X3+B2         (B2) = LOC(CXI1 )
*                              (I=1)
 ONE      SA2    B2            (X2) = REAL(CXI1 )
          SA3    B2-B7         (X3) = IMAG(CXI1 )
          BX0    X2
          BX1    X3
          AX2    59
          AX3    59
          BX4    X2-X0         (X4) = ABS(REAL(CXI1 ))
          BX5    X3-X1         (X5) = ABS(IMAG(CXI1 ))
*
          RX0    X6+X4
          NX0    X0
          RX1    X0+X5         (X6) = (X6) + (X5) + (X4)
          NX6    X1
*
          ZR     B1,OUT        IF I .EQ. N , GO TO OUT
*
*                              (I = I+1)
 LOOP     SA2    A2+B3         (X2) = REAL(CXII )
          SA3    A2-B7         (X3) = IMAG(CXII )
          BX0    X2
          BX1    X3
          AX2    59
          AX3    59
          BX4    X2-X0         (X4) = ABS(REAL(CXII ))
          BX5    X3-X1         (X5) = ABS(IMAG(CXII ))
*
          RX0    X6+X4         (X6) = (X6) + (X5) + (X4)
          NX0    X0
          RX1    X0+X5
          SB1    B1+B7         COUNT TERM
          NX6    X1
*
          NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
 OUT      OUTFTN SCASUM        RETURN
          END
*DECK,SSCAL
          IDENT  SSCAL
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL SSCAL(N,SA,SX,INCX)
*
*         SA*SXII   REPLACES SXII   FOR I=1,N
*
*         SXII  = SX(1 + (I-1)*INCX)  IF INCX .GE. 0
*               = SX(1 + (I-N)*INCX)  IF INCX .LT. 0
*
*         SX( )                     SINGLE PRECISION
*         N,INCX                    INTEGER TYPE
*         SA                        SINGLE PRECISION
*
*         ROUNDED ARITHMETIC INSTRUCTIONS ARE USED
*
*         WRITTEN BY  DAVID R. KINCAID AND ELIZABETH WILLIAMS
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  SSCAL
          VFD    42/5HSSCAL,18/4
*
 SSCAL    DATA   0             ENTRY/EXIT
          INFTN  SSCAL,4
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          SB1    X1+B7         (B1) = N-1
          SA2    B2            (X2) = SA
*
          SA4    B4            (X4) = INCX
          NG     B1,OUT        IF N .LE. 0 , GO TO OUT
          SX1    -B1           (X1) = -(N-1)
          SB4    X4            (B4) = INCX
*
          GT     B4,ONE        IF INCX .GT. 0 , GO TO ONE
          DX4    X1*X4         LOC(SXI1 ) = LOC(SX) - (N-1)*INCX
          SB3    X4+B3         (B3) = LOC(SXI1 )
*
*                              (I = 1)
 ONE      SA3    B3            (X3) = SXI1
          FX6    X2*X3         (X6) = SA*SXI1
          SA6    B3
*
          ZR     B1,OUT        IF I .EQ. N , GO TO OUT
*
*                              (I = I+1)
 LOOP     SA3    A3+B4         (X3) = SXII
          FX6    X2*X3         (X6) = SA*SXII
          SB1    B1+B7         I = I+1
          SA6    A3            SXII  = (X6)
*
          NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
 OUT      OUTFTN SSCAL         RETURN
          END
*DECK,DSCAL
          IDENT  DSCAL
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL DSCAL(N,DA,DX,INCX)
*
*         DA*DXII   REPLACES DXII   FOR I=1,N
*
*         DXII  = DX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*               = DX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         DX( )                     DOUBLE PRECISION
*         N,INCX                    INTEGER TYPE
*         DA                        DOUBLE PRECISION
*
*         WRITTEN BY  DAVID R. KINCAID AND ELIZABETH WILLIAMS
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  DSCAL
          VFD    42/5HDSCAL,18/4
*
 DSCAL    DATA   0             ENTRY/EXIT
          INFTN  DSCAL,4
          SA3    B1            (X3) = N
          SB7    -1            (B7) = -1
          SB1    X3+B7         (B1) = N-1
          SA1    B2            (X1,X2) = DA
          SA2    B2-B7
*
          SA4    B4            (X4) = INCX
          NG     B1,OUT        IF N .LE. 0 , GO TO OUT
          LX4    1             INCX = 2*INCX
          SX3    -B1           (X3) = -(N-1)
          SB4    X4            (B4) = INCX
*
          GT     B4,ONE        IF INCX .GT. 0 , GO TO ONE
          DX4    X3*X4         LOC(DXI1 ) = LOC(DX) - (N-1)*INCX
          SB3    X4+B3
*
*                              (I = 1)
 ONE      SA3    B3            (X3,X4) = DXI1
          SA4    B3-B7
*
          FX5    X2*X3         (X6,X7) = DA*DXI1
          FX0    X1*X4
          FX5    X0+X5
          FX4    X1*X3
          DX0    X1*X3
          FX5    X0+X5
          FX6    X4+X5
          DX7    X4+X5
*
          SA6    A3            DXI1  = (X6,X7)
          SA7    A4
*
          ZR     B1,OUT        IF I .EQ. N , GO TO OUT
*
*                              (I = I+1)
 LOOP     SA3    A3+B4         (X3,X4) = DXII
          SA4    A3-B7
*
*
          FX5    X2*X3         (X6,X7) = DA*DXII
          FX0    X1*X4
          FX5    X0+X5
          FX4    X1*X3
          DX0    X1*X3
          FX5    X0+X5
          FX6    X4+X5
          DX7    X4+X5
*
          SB1    B1+B7         I = I+1
*
          SA6    A3            DXII  = (X6,X7)
          SA7    A4
*
          NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
 OUT      OUTFTN DSCAL         RETURN
          END
*DECK,CSCAL
          IDENT  CSCAL
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL CSCAL(N,CA,CX,INCX)
*
*         CA*CXII   REPLACES CXII   FOR I=1,N
*
*         CXII  = CX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*               = CX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         CX( )                     COMPLEX TYPE
*         N,INCX                    INTEGER TYPE
*         CA                        COMPLEX TYPE
*
*         ROUNDED ARITHMETIC INSTRUCTIONS ARE USED
*
*         WRITTEN BY  DAVID R. KINCAID AND ELIZABETH WILLIAMS
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  CSCAL
          VFD    42/5HCSCAL,18/4
*
 CSCAL    DATA   0             ENTRY/EXIT
          INFTN  CSCAL,4
          SA3    B1            (X3) = N
          SB7    -1            (B7) = -1
          SB1    X3+B7         (B1) = N-1
          SA1    B2            (X1) = REAL(CA)
          SA2    B2-B7         (X2) = IMAG(CA)
*
          NG     B1,OUT        IF N .LE. 0 , GO TO OUT
          SA4    B4            (X4) = INCX
          LX4    1             INCX = 2*INCX
          SX3    -B1           (X3) = -(N-1)
          SB4    X4            (B4) = INCX
*
          GT     B4,ONE        IF INCX .GT. 0 , GO TO ONE
          DX4    X3*X4         LOC(CXI1 ) = LOC(CX) - (N-1)*INCX
          SB3    X4+B3         (B3) = LOC(CXI1 )
*
*                              (I = 1)
 ONE      SA3    B3            (X3) = REAL(CXI1 )
          SA4    B3-B7         (X4) = IMAG(CXI1 )
*
*                              (X6,X7) = CA*CXI1
          RX6    X1*X3         (X6) = REAL(CA)*REAL(CXI1 )
          RX5    X2*X4         (X5) = IMAG(CA)*IMAG(CXI1 )
          RX0    X6-X5         (X0) = REAL(CA*CXI1 )
          NX6    X0
*
          RX7    X1*X4         (X7) = REAL(CA)*IMAG(CXI1 )
          RX5    X2*X3         (X5) = IMAG(CA)*REAL(CXI1 )
          RX0    X7+X5         (X0) = IMAG(CA*CXI1 )
          NX7    X0
*
          SA6    A3            REAL(CXI1 ) = (X6)
          SA7    A4            IMAG(CXI1 ) = (X7)
*
          ZR     B1,OUT        IF I .EQ. N , GO TO OUT
*
*                              (I = I+1)
 LOOP     SA3    A3+B4         (X3) = REAL(CXII )
          SA4    A3-B7         (X4) = IMAG(CXII )
*
          RX6    X1*X3         (X6) = REAL(CA)*REAL(CXII )
          RX5    X2*X4         (X5) = IMAG(CA)*IMAG(CXII )
          RX0    X6-X5         (X0) = REAL(CA*CXII )
          NX6    X0
*
          RX7    X1*X4         (X7) = REAL(CA)*IMAG(CXII )
          RX5    X2*X3         (X5) = IMAG(CA)*REAL(CXII )
          RX0    X7+X5         (X0) = IMAG(CA*CXII )
          SB1    B1+B7         I = I+1
          NX7    X0
*
          SA6    A3            REAL(CXII ) = (X6)
          SA7    A4            IMAG(CXII ) = (X7)
*
          NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
 OUT      OUTFTN CSCAL         RETURN
          END
*DECK,CSSCAL
          IDENT  CSSCAL
*
***       USE WITH FORTRAN STATEMENT
*
*         CALL CSSCAL(N,SA,CX,INCX)
*
*         SA*CXII   REPLACES CXII   FOR I=1,N
*
*         CXII  = CX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*               = CX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         CX( )                     COMPLEX TYPE
*         N,INCX                    INTEGER TYPE
*         SA                        SINGLE PRECISION
*
*         ROUNDED ARITHMETIC INSTRUCTIONS ARE USED
*
*         WRITTEN BY  DAVID R. KINCAID AND ELIZABETH WILLIAMS
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  CSSCAL
          VFD    42/6HCSSCAL,18/4
*
 CSSCAL   DATA   0             ENTRY/EXIT
          INFTN  CSSCAL,4
          SA3    B1            (X3) = N
          SB7    -1            (B7) = -1
          SB1    X3+B7         (B1) = N-1
          SA2    B2            (X2) = SA
*
          SA4    B4            (X4) = INCX
          NG     B1,OUT        IF N .LE. 0 , GO TO OUT
          LX4    1             INCX = 2*INCX
          SX3    -B1           (X3) = -(N-1)
          SB4    X4            (B4) = INCX
*
          GT     B4,ONE        IF INCX .GT. 0 , GO TO ONE
          DX4    X3*X4         LOC(CXI1 ) = LOC(CX) - (N-1)*INCX
          SB3    X4+B3         (B3) = LOC(CXI1 )
*
*                              (I = 1)
 ONE      SA3    B3            (X3) = REAL(CXI1 )
          SA4    B3-B7         (X4) = IMAG(CXI1 )
*
*                              (X6,X7) = SA*CXI1
          RX6    X2*X3         (X6) = SA*REAL(CXI1 )
          RX7    X2*X4         (X7) = SA*IMAG(CXI1 )
*
          SA6    A3            REAL(CXI1 ) = (X6)
          SA7    A4            IMAG(CXI1 ) = (X7)
*
          ZR     B1,OUT        IF I .EQ. N , GO TO OUT
*
*                              (I = I+1)
 LOOP     SA3    A3+B4         (X3) = REAL(CXII )
          SB1    B1+B7         I = I+1
          SA4    A3-B7         (X4) = IMAG(CXII )
*
*                              (X6,X7) = SA*CXII
          RX6    X2*X3         (X6) = SA*REAL(CXII )
          RX7    X2*X4         (X7) = SA*IMAG(CXII )
*
          SA6    A3            REAL(CXII ) = (X6)
          SA7    A4            IMAG(CXII ) = (X7)
*
          NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
 OUT      OUTFTN CSSCAL        RETURN
          END
*DECK,ISAMAX
          IDENT  ISAMAX
*
***       INTEGER FUNCTION ISAMAX(N,SX,INCX)
*
*         FIND AN INDEX  I(MAX)  CORRESPONDING TO THE MAXIMUM ABSOLUTE V
*         COMPONENTS  SXII   OF THE VECTOR SX.
*
*         SXII  = SX(1 + (I-1)*INCX)  IF INCX .GE. 0
*               = SX(1 + (I-N)*INCX)  IF INCX .LT. 0
*
*         SX( )                     SINGLE PRECISION
*         N,INCX                    INTEGER TYPE
*         RESULT ISAMAX             INTEGER TYPE
*
*         WRITTEN BY  DAVID R. KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  ISAMAX
          VFD    42/6HISAMAX,18/3
*
 ISAMAX   DATA   0             ENTRY/EXIT
          INFTN  ISAMAX,3
          MX6    0             (X6)=ISAMAX=0
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          SB4    X1            (B4) = N
          SB1    X1+B7         (B1) = N-1
          NG     B1,OUT        IF(N .LE. 0) GO TO OUT
*
          SX6    -B7           (X6) = 1             (ISAMAX)
          LE     B1,OUT        IF N .LE. 1 , GO TO OUT
*
          SA3    B3            (X3) = INCX
          SX1    -B1           (X1) = -(N-1)
          SB3    X3            (B3) = INCX
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3         LOC(XI1 ) = LOC(SX) - (N-1)*INCX
          SB2    X3+B2         (B2) = LOC(SXI1 )
*
*                              (I = 1)
 ONE      SA2    B2            (X2) = SXI1
          BX3    X2
          AX2    59
          BX5    X2-X3         (X5) = ABS(SXI1 )    (SAMAX)
*
*
*                              (I=I+1)
 LOOP     SA2    A2+B3         (X2) = SXII
          BX3    X2
          AX2    59
          BX2    X2-X3         (X2) = ABS(SXII )
          SB1    B1+B7         COUNT TERM
*
          NX5    X5
          FX0    X5-X2
          PL     X0,TEST       IF ABS(SXII ) .LE. SAMAX , GO TO TEST
*
          BX5    X2            (X5) = ABS(SXII )    (SAMAX)
          SX6    B4-B1         (X6) = I             (ISAMAX)
*
 TEST     NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
 OUT      OUTFTN ISAMAX        RETURN
          END
*DECK,IDAMAX
          IDENT  IDAMAX
*
***       INTEGER FUNCTION IDAMAX(N,DX,INCX)
*
*         FIND AN INDEX  I(MAX)  CORRESPONDING TO THE MAXIMUM ABSOLUTE V
*         COMPONENTS  DXII   OF THE VECTOR   DX
*
*         DXII  = DX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*               = DX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         DX( )                     DOUBLE PRECISION
*         N,INCX                    INTEGER TYPE
*         RESULT IDAMAX             INTEGER TYPE
*
*         WRITTEN BY  DAVID R. KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  IDAMAX
          VFD    42/6HIDAMAX,18/3
*
 IDAMAX   DATA   0             ENTRY/EXIT
          INFTN  IDAMAX,3
          MX6    0             (X6)=IDAMAX=0
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          SB4    X1            (B4) = N
          SB1    X1+B7         (B1) = N-1
          NG     B1,OUT        IF(N .LE. 0) GO TO OUT
*
          SX6    -B7           (X6) = 1
          LE     B1,OUT        IF N .LE. 1 , GO TO OUT
*
          SA3    B3            (X3) = INCX
          SX1    -B1           (X1) = -(N-1)
          LX3    1             INCX = 2*INCX
          SB3    X3            (B3) = INCX
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3         LOC(DXI1 ) = LOC(DX) - (N-1)*INCX
          SB2    X3+B2         (B2) = LOC(DXI1 )
*
*                              (I=1)
 ONE      SA2    B2
          SA3    B2-B7         (X2,X3) = DXI1
          BX0    X2
          AX0    59
          BX4    X0-X2
          BX5    X0-X3         (X4,X5) = DABS(DXI1 )
*
*                              (I=I+1)
 LOOP     SA2    A2+B3
          SA3    A3+B3         (X2,X3) = DXII
          BX0    X2
          AX0    59
          BX2    X0-X2
          BX3    X0-X3         (X2,X3) = DABS(DXII )
          SB1    B1+B7         COUNT TERM
*
          FX1    X4-X2         IF DABS(DXII ) .LE. DAMAX , GO TO TEST
          FX5    X5-X3
          DX4    X4-X2
          NX1    X1
          FX4    X4+X5
          NX5    X4
          FX4    X1+X5
          PL     X4,TEST
*
          SX6    B4-B1         (X6) = I                (IDAMAX)
          BX4    X2            (X4,X5) = DABX(DXII )   (DAMAX)
          BX5    X3
*
 TEST     NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
 OUT      OUTFTN IDAMAX        RETURN
          END
*DECK,ICAMAX
          IDENT  ICAMAX
*
***       INTEGER FUNCTION ICAMAX(N,CX,INCX)
*
*         FIND AN INDEX  I(MAX)  CORRESPONDING TO THE MAXIMUM SUM OF THE
*         ABSOLUTE VALUE OF THE REAL PART AND THE ABSOLUTE VALUE OF THE
*         IMAGINARY PART OF THE COMPONENTS  CXII   OF THE VECTOR CX
*
*         CXII  = CX(1 + (I-1)*2*INCX)  IF INCX .GE. 0
*               = CX(1 + (I-N)*2*INCX)  IF INCX .LT. 0
*
*         CX( )                     COMPLEX TYPE
*         N,INCX                    INTEGER TYPE
*         RESULT ICAMAX             INTEGER TYPE
*
*         WRITTEN BY  DAVID R. KINCAID
*                     CENTER FOR NUMERICAL ANALYSIS/COMPUTATION CENTER
*                     THE UNIVERSITY OF TEXAS AT AUSTIN
***       1 JUNE 77
*
          ENTRY  ICAMAX
          VFD    42/6HICAMAX,18/3
*
 ICAMAX   DATA   0             ENTRY/EXIT
          INFTN  ICAMAX,3
          MX6    0             (X6)=ICAMAX=0
          SA1    B1            (X1) = N
          SB7    -1            (B7) = -1
          SB4    X1            (B4) = N
          SB1    X1+B7         (B1) = N-1
          NG     B1,OUT        IF(N .LE. 0) GO TO OUT
*
          SX6    -B7           (X6) = 1
          LE     B1,OUT        IF N .LE. 1 , GO TO OUT
*
          SA3    B3            (X3) = INCX
          SX1    -B1           (X1) = -(N-1)
          LX3    1             (X3) = 2*INCX
          SB3    X3            (B3) = INCX
*
          GT     B3,ONE        IF INCX .GT. 0 , GO TO ONE
          DX3    X1*X3         LOC(CXI1 ) = LOC(CX) - (N-1)*INCX
          SB2    X3+B2         (B2) = LOC(CXI1 )
*
*                              (I = 1)
 ONE      SA2    B2            (X2) = REAL(CXI1 )
          BX3    X2
          AX2    59
          BX5    X2-X3         (X5) = ABS(REAL(CXI1 ))
          SA3    B2-B7         (X3) = IMAG(CXI1 )
          BX2    X3
          AX3    59
          BX4    X3-X2         (X4) = ABS(IMAG(CXI1 )
*
          RX5    X4+X5
          NX5    X5            (X5) = (X4) + (X5)    (AMAX)
*
*                              (I = I+1)
 LOOP     SA2    A2+B3         (X2) = REAL(CXII )
          BX3    X2
          AX2    59
          BX2    X2-X3         (X2) = ABS(REAL(CXII ))
          SA3    A2-B7         (X3) = IMAG(CXII )
          BX7    X3
          AX3    59
          BX3    X3-X7         (X3) = ABS(IMAG(CXII )
*
          RX2    X2+X3
          SB1    B1+B7         COUNT TERM
          NX2    X2            (X2) = (X2) + (X3)
*
          FX0    X5-X2
          PL     X0,TEST       IF  ABS(REAL(CXII )) + ABS(IMAG(CXII )) .
*
          BX5    X2            (X5) = ABS(REAL(CXII ))    (AMAX)
          SX6    B4-B1         (X6) = I    (ICAMAX)
*
 TEST     NZ     B1,LOOP       IF I .NE. N , GO TO LOOP
*
 OUT      OUTFTN ICAMAX        RETURN
          END
          AXR$
