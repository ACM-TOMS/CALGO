      SUBROUTINE ZZBSTR ( STRING, START, END, PTR, NEWSTR, MAX, IGNORE )

C## A R G U M E N T S:
                       CHARACTER * (*) STRING,  NEWSTR
                       INTEGER         PTR,  START,  END
                       LOGICAL         MAX, IGNORE

C## S T A T U S:
C               SINGLE/DOUBLE CONVERSION:        NOT REQUIRED.
C               SYSTEM  DEPENDENCE:                      NONE.
C
C>RCS $HEADER: BSTR.F,V 1.10 91/11/20 10:52:38 BUCKLEY EXP $
C>RCS $LOG:     BSTR.F,V $
C>RCS REVISION 1.10  91/11/20  10:52:38  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 1.9  89/06/30  13:38:13  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3.1.1  89/05/20  16:42:31  BUCKLEY
C>RCS TEMP. TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.3  89/05/18  12:20:17  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:47:06  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C>RCS REVISION 1.1  89/01/17  16:50:06  BUCKLEY
C>RCS INITIAL REVISION
C>RCS
C
C## D E S C R I P T I O N:
C
C     THIS IS FOR BUILDING A STRING IN "STRING".  IT ADDS CHARAC-
C     TERS FROM NEWSTR.  IT ASSUMES THAT ON ENTRY, STRING CONTAINS
C     PTR CHARACTERS; PTR IS UPDATED BEFORE EXIT.
C
C     IGNORE CONTROLS THE PROCESSING OF TRAILING BLANKS IN NEWSTR.
C     IF IGNORE IS TRUE, THEY ARE STRIPPED OFF, I.E. IGNORED;
C     OTHERWISE THEY ARE COPIED TO STRING JUST LIKE ANY OTHER
C     CHARACTER.
C
C     START AND END CONTROL THE PLACEMENT OF NEWSTR IN STRING.
C     IF START > PTR ON ENTRY, NEWSTR IS APPENDED ONTO STRING,
C     STARTING IN POSITION PTR+1.  IF START <= PTR ON ENTRY, THEN
C     NEWSTR IS OVERWRITTEN INTO STRING IN POSITIONS START TO END.
C     IF THERE IS NOT ENOUGH ROOM, THEN THE TAIL OF STRING IS RIGHT
C     SHIFTED JUST ENOUGH TO MAKE ROOM. WHEN END < START, THE NEW
C     STRING IS INSERTED JUST BEFORE POSITION START.
C
C     MAX CONCERNS THE ACTION TAKEN WHEN STRING OVERFLOWS.  IF MAX IS
C     FALSE ON ENTRY, ANY EXCESS CHARACTERS ARE JUST IGNORED AND AS
C     MANY AS POSSIBLE ARE PUT INTO STRING. PTR IS SET TO THE LENGTH
C     OF STRING.  IF MAX IS TRUE, THEN, WHEN THERE ARE TOO MANY
C     CHARACTERS IN NEWSTR TO FIT ON THE END OF STRING, THERE IS NO
C     ACTION TAKEN, AND BOTH STRING AND PTR ARE UNCHANGED. IN EITHER
C     CASE, ON EXIT, MAX IS TRUE IF STRING DID NOT HAVE ENOUGH ROOM
C     TO ACCEPT ALL OF THE REQUIRED CHARACTERS (INCLUDING THE
C     BLANKS IF APPROPRIATE) FROM NEWSTR.
C
C## E N T R Y   P O I N T S: THE NATURAL ENTRY ZZBSTR
C## S U B R O U T I N E S:   ZZLENG   NONBLANK LENGTH OF A STRING.
C                            ZZSHFT   SHIFT A STRING.
C                            LEN   ...INTRINSIC
C## P A R A M E T E R S:     NONE ARE DEFINED.
C## L O C A L   D E C L:
                             INTEGER  LENS, LENN, ZZLENG, K, ROOM

C## S A V E:                 NONE SELECTED.
C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:                 NONE ARE SET.
C##                                                E X E C U T I O N
C##                                                E X E C U T I O N

      LENS = LEN ( STRING )

      IF ( IGNORE ) THEN
         LENN = ZZLENG ( NEWSTR )
      ELSE
         LENN =   LEN  ( NEWSTR )
      ENDIF

      IF ( START .LE. PTR ) THEN
         IF ( START .GT. END ) THEN
            ROOM = 0
            K    = -LENN
         ELSE
            ROOM = END - START + 1
            K    = ROOM - LENN
         ENDIF

         IF ( K .GT. 0 ) THEN
            STRING(START:END) = NEWSTR(1:LENN)
            MAX = .FALSE.
         ELSE
            ROOM = PTR - K
            IF ( ROOM .GT. LENS ) THEN
               IF ( MAX ) THEN
                  GOTO 90000
               ELSE
                  ROOM = LENS
                  MAX  = .TRUE.
               ENDIF
            ELSE
               MAX = .FALSE.
            ENDIF

            IF ( END .GT. START ) THEN
               CALL ZZSHFT ( STRING, END+1, END-K+1, LENS )
            ELSE
               CALL ZZSHFT(STRING,START,START+LENN,LENS)
            ENDIF
            STRING(START:START+LENN -1) = NEWSTR(1:LENN)
            PTR = ROOM
         ENDIF

      E L S E

         K = PTR + LENN
         IF ( K .LE. LENS ) THEN
            STRING(PTR+1:K) = NEWSTR(1:LENN)
            MAX = .FALSE.
            PTR =  K
         ELSE
            IF ( MAX ) THEN
               GOTO 90000
            ELSE IF ( PTR .LT. LENS ) THEN
               STRING ( PTR+1 : LENS ) = NEWSTR ( 1:LENN )
               PTR = LENS
               MAX = .TRUE.
            ENDIF
         ENDIF
      ENDIF

      PTR = ZZLENG ( STRING )

C## E X I T
90000       RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZBSTR.
                    END