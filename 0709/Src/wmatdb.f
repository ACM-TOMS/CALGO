      SUBROUTINE ZZWMAT (SPACE, NAME, A,  X,  B,  Y,
     -                   NRX, R1, R2, C1, C2,
     -                   F, W, D, S, LINE, UNIT )

C## A R G U M E N T S:
                        CHARACTER *(*) NAME
                        CHARACTER *1   F
                        INTEGER  NRX, R1, R2, C1, C2
                        INTEGER  W, D, UNIT, LINE, S, SPACE

                        DOUBLE PRECISION X(NRX,*), Y(NRX,*), A, B
C!!!!                   REAL             X(NRX,*), Y(NRX,*), A, B

C## S T A T U S:
C               SINGLE/DOUBLE CONVERSION: NEEDED (SEE CONVRT).
C
C               IGNORE LINES BEGINNING WITH  "C!!!!" .
C
C               THIS VERSION IS IN   D O U B L E   PRECISION.
C!!!!           THIS VERSION IS IN   S I N G L E   PRECISION.
C
C               SYSTEM  DEPENDENCE:                      NONE.

C>RCS $HEADER: WMAT.F,V 1.3 91/12/16 11:57:27 BUCKLEY EXP $
C>RCS $LOG:     WMAT.F,V $
C>RCSREVISION 1.3  91/12/16  11:57:27  BUCKLEY
C>RCSFIX FOR TOMS; ADD SPACE.
C>RCS
C>RCSREVISION 1.1  90/07/31  12:57:12  BUCKLEY
C>RCSINITIAL REVISION
C>RCS

C## D E S C R I P T I O N:

C       PRINT OUT A SUBMATRIX OF A MATRIX WITH GIVEN FORMAT.

C       PRINT A TITLE FOR THE MATRIX: NAME

C       EACH OUTPUT LINE IS PRECEDED BY  SPACE  BLANK CHARACTERS.

C       PRINT THE VALUES FROM THE MATRIX A*X + B*Y
C         A AND B ARE SCALARS; X AND Y ARE VECTORS OR MATRICES.
C         THERE ARE ASSUMED TO BE NRX ROWS IN EACH MATRIX.
C           (USE 1 TO PRINT COLUMN VECTORS.)

C       PRINT FROM ROW    R1 TO ROW    R2, INCLUSIVE.
C       PRINT FROM COLUMN C1 TO COLUMN C2, INCLUSIVE.

C       PRINT EACH ENTRY IN FORMAT  FW.D
C              ..F IS A CHARACTER 'F', 'G', 'E' OR 'D'
C              ..W AND D ARE INTEGERS
C              ..S IS THE NUMBER OF SPACES BETWEEN EACH ENTRY.
C              ..LINE IS THE NUMBER OF CHARACTERS PER LINE.
C                  (IF LINE=0, THEN LINE IS REPLACED WITH 80)

C       ALL OUTPUT IS ON THE GIVEN UNIT.

C## E N T R Y   P O I N T S: THE NATURAL ENTRY ZZWMAT.
C## S U B R O U T I N E S:   NONE ARE CALLED.
C## P A R A M E T E R S:     NONE ARE DEFINED.
C## L O C A L   D E C L:

        INTEGER  I, J

        CHARACTER *100 FORM

C## S A V E:                 NONE SELECTED.
C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:

        DATA FORM / '(??X,??(? ??.??, ?? X), ? ??.??)' /
C        E.G.        (02X,05(F 09.04, 03 X), F 09.04)
C                    123456789012345678901234567890
C                             1         2         3

C##                                                E X E C U T I O N
C##                                                E X E C U T I O N

        IF ( LINE .EQ. 0 ) THEN
            LINE = 80
        ENDIF

        IF ( SPACE .EQ. 0 ) THEN
            WRITE ( FORM(02:05), 99998 ) ' '
        ELSE
            WRITE ( FORM(02:03), 99999 ) SPACE
        ENDIF
        WRITE ( FORM(06:07), 99999 ) 1+(LINE-SPACE-W)/(W+S)
        WRITE ( FORM(09:09), 99997 ) F
        WRITE ( FORM(11:12), 99999 ) W
        WRITE ( FORM(14:15), 99999 ) D
        WRITE ( FORM(18:19), 99999 ) S
        WRITE ( FORM(25:25), 99997 ) F
        WRITE ( FORM(27:28), 99999 ) W
        WRITE ( FORM(30:31), 99999 ) D

        IF ( SPACE .NE. 0 ) THEN
            WRITE ( UNIT, 99998 ) (' ',I=1,SPACE),NAME
        ELSE
            WRITE ( UNIT, 99998 ) NAME
        ENDIF

        DO 1000 I = R1, R2

            IF ( SPACE .NE. 0 ) THEN
                WRITE ( UNIT, FORM )  SPACE,(A*X(I,J)+B*Y(I,J),J=C1,C2)
            ELSE
                WRITE ( UNIT, FORM )  ( A*X(I,J)+B*Y(I,J), J=C1, C2)
            ENDIF

 1000   CONTINUE

C## E X I T
90000      RETURN

C## F O R M A T S:  SEE CHARACTER STRING FORM.

99997 FORMAT (A1)
99998 FORMAT (A)
99999 FORMAT (I2)

C                 E N D      OF WMAT.
      END
