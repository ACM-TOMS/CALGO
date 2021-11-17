      SUBROUTINE ZZSECS(SECS)

C## A R G U M E N T S:
                      DOUBLE PRECISION  SECS
C!!!!                 REAL              SECS

C## S T A T U S:
C      SINGLE/DOUBLE CONVERSION:      NEEDED (SEE CONVRT).
C
C      IGNORE LINES BEGINNING WITH  "C!!!!" .
C
C      THIS VERSION IS IN   D O U B L E   PRECISION.
C!!!!  THIS VERSION IS IN   S I N G L E   PRECISION.
C
C      SYSTEM  DEPENDENCE:   SYSTEM ROUTINE FOR CPU USAGE.

C             THIS VERSION IS FOR  SUN4
C
C>RCS $HEADER: SECS.GL,V 2.1 91/11/22 11:45:25 BUCKLEY EXP $
C>RCS $LOG:     SECS.GL,V $
C>RCS REVISION 2.1  91/11/22  11:45:25  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 2.0  90/07/06  10:48:10  BUCKLEY
C>RCS COMMON VERSION FOR TOMS AND MT
C>RCS
C>RCS REVISION 1.9  89/06/30  13:30:19  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3.1.1  89/05/20  13:46:31  BUCKLEY
C>RCS TEMP. TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.3  89/05/18  12:13:31  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:35:12  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C>RCS REVISION 1.1  89/01/17  16:34:33  BUCKLEY
C>RCS INITIAL REVISION
C>RCS
C
C## D E S C R I P T I O N:

C     THIS SUBROUTINE SHOULD OBTAIN (FROM THE OPERATING SYSTEM) THE
C     AMOUNT OF CPU TIME USED BY THE CALLING PROGRAM SINCE THE
C     EXECUTION BEGAN. IF DESIRABLE, "SECS" CAN ALSO BE CONVERTED
C     TO DOUBLE PRECISION (SEE CONVRT). HOWEVER, THE ROUTINE ACTUALLY
C     WORKS TOTALLY AS A SINGLE PRECISION ROUTINE, EXCEPT THAT THE
C     VALUE WHICH IS PASSED BACK MAY BE IN EITHER PRECISION AS
C     APPROPRIATE.

C     TIME IS MEASURED FROM THE FIRST CALL TO ZZSECS.  THUS
C     ON THE FIRST CALL TO ZZSECS, A TIME OF  0.0  SECONDS IS ALWAYS
C     RETURNED.

C## E N T R Y   P O I N T S: THE NATURAL ENTRY ZZSECS.
C## S U B R O U T I N E S:   A SYSTEM CLOCK.

C## P A R A M E T E R S:


      REAL        ZERO
      PARAMETER ( ZERO = 0.0E0 )

C## L O C A L   D E C L:

      LOGICAL  FIRST


      REAL    ETIME, DUMMY(2)
      REAL    STTIME,  SEC


C## S A V E:

      SAVE   FIRST, STTIME

C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.

C## D A T A:

      DATA   FIRST/.TRUE./

C##                                            E X E C U T I O N
C##                                            E X E C U T I O N

      IF ( FIRST ) THEN

         FIRST  = .FALSE.

         STTIME = ETIME(DUMMY)

         SEC   = ZERO

      ELSE

         SEC = ETIME(DUMMY) - STTIME

      ENDIF

      SECS = DBLE(SEC)
C!!!! SECS =      SEC

C## E X I T
90000      RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZSECS.
                    END
