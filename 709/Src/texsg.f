      SUBROUTINE ZZTEX (UNIT, NAME, PNO, FVALS, ITERS, FC, DF, FXSTAR,
     -                  ACC, GNORM, GMIN, DISTX, CRIT, ERR, FTIME,MTIME)

C## A R G U M E N T S:
                        INTEGER  UNIT, PNO, FVALS, ITERS, CRIT, ERR,ACC
                        REAL             FC, FXSTAR, GNORM, GMIN, DISTX
C!!!!                   DOUBLE PRECISION FC, FXSTAR, GNORM, GMIN, DISTX
                        REAL             DF, FTIME, MTIME
C!!!!                   DOUBLE PRECISION DF, FTIME, MTIME
                        CHARACTER *(*) NAME

C## S T A T U S:
C               IGNORE LINES BEGINNING WITH  "C!!!!" .
C
C               SYSTEM  DEPENDENCE:                      NONE.
C               SINGLE/DOUBLE CONVERSION: NEEDED (SEE CONVRT).
C
C               THIS VERSION IS IN   S I N G L E   PRECISION.
C!!!!           THIS VERSION IS IN   D O U B L E   PRECISION.
C
C>RCS $HEADER: TEX.F,V 2.1 91/11/20 10:53:11 BUCKLEY EXP $
C>RCS $LOG:     TEX.F,V $
C>RCS REVISION 2.1  91/11/20  10:53:11  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 2.0  90/07/31  11:37:38  BUCKLEY
C>RCS MINOR FORMAT FIX.
C>RCS
C>RCS REVISION 1.9.1.1  89/07/02  14:36:13  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.9  89/06/30  13:39:54  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3.1.1  89/05/20  16:43:02  BUCKLEY
C>RCS TEMP. TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.3  89/05/18  12:20:55  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:48:23  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C
C## D E S C R I P T I O N:

C     THIS ROUTINE PRODUCES THE "OUTPUT" FILE OF COMPUTED SOLUTIONS
C     FOR EACH TEST FUNCTION.

C## E N T R Y   P O I N T S: THE NATURAL ENTRY ZZTEX.
C## S U B R O U T I N E S:   NONE ARE CALLED.
C## P A R A M E T E R S:     NONE ARE DEFINED.


      LOGICAL     T,          F
      PARAMETER ( T = .TRUE., F = .FALSE. )

      CHARACTER*(*) TRUE,          QT,       FALSE,           QF
      PARAMETER (   TRUE = 'TRUE', QT = 'T', FALSE = 'FALSE', QF = 'F' )

      INTEGER     ITRUE,     IFALSE
      PARAMETER ( ITRUE = 1, IFALSE = 0 )

      REAL              RTRUE,        RFALSE
C!!!! DOUBLE PRECISION  RTRUE,        RFALSE
      PARAMETER      (  RTRUE = 1.D0, RFALSE = 0.D0 )

      INTEGER     CTOUPP,     CTOLOW,     CTOCAP
      PARAMETER ( CTOUPP = 1, CTOLOW = 2, CTOCAP = 3 )

      CHARACTER*(*) BLANK,        QUOTE,        HASH
      PARAMETER (   BLANK  = ' ', QUOTE  = '"', HASH   = '#' )

      CHARACTER*(*) PERIOD,       COMMA,        SEMICN
      PARAMETER (   PERIOD = '.', COMMA  = ',', SEMICN = ';' )

      CHARACTER*(*) COLON,        DASH,         EQUALS
      PARAMETER (   COLON  = ':', DASH   = '-', EQUALS = '=' )

      CHARACTER*(*) OBRACE,       CBRACE,       UNDERS
      PARAMETER (   OBRACE = '{', CBRACE = '}', UNDERS = '_' )

      CHARACTER*(*) PLUS,         MINUS,        EXCLAM
      PARAMETER (   PLUS   = '+', MINUS  = '-', EXCLAM = '!' )

      CHARACTER*(*) GTHAN,        LTHAN,        QUESMK
      PARAMETER (   GTHAN  = '>', LTHAN  = '<', QUESMK = '?' )

      CHARACTER*(*) SLASH,        BSLASH,       PERCNT
      PARAMETER (   SLASH  = '/', BSLASH = '\\',PERCNT = '%' )

      CHARACTER*(*) CARAT,        ATSIGN,       TILDE
      PARAMETER (   CARAT  = '^', ATSIGN = '@', TILDE = '~' )
      REAL              ZERO,       ONE,       TWO,       THREE
C!!!! DOUBLE PRECISION  ZERO,       ONE,       TWO,       THREE
      PARAMETER (       ZERO = 0D0, ONE = 1D0, TWO = 2D0, THREE = 3D0)

      REAL              FOUR,       FIVE,      SIX,       SEVEN
C!!!! DOUBLE PRECISION  FOUR,       FIVE,      SIX,       SEVEN
      PARAMETER (       FOUR = 4D0, FIVE = 5D0,SIX = 6D0, SEVEN = 7D0)

      REAL              EIGHT,         NINE,          TEN
C!!!! DOUBLE PRECISION  EIGHT,         NINE,          TEN
      PARAMETER (       EIGHT = 8D0,   NINE = 9D0,    TEN = 10D0     )

        CHARACTER *(*) COMPTD, TEX, FP, NOX, PATH
        PARAMETER ( COMPTD = BSLASH//'COMPUTED'//PERCNT,
     -              TEX    = '.TEX',
     -              NOX    = '   ',
     -              PATH   = 'COMPUTED/',
     -              FP     = BSLASH//'FP'
     -  )
        INTEGER     LINELN,       FNAMLN
        PARAMETER ( LINELN = 100, FNAMLN = 80 )

C## L O C A L   D E C L:
                         CHARACTER *(LINELN) OUT(3)
                         CHARACTER *(FNAMLN) FNAME, LSTNAM
                         INTEGER             POS, LASTPR, POS1, POS2
                         INTEGER             ZZLENG, ZZLFTI
                         LOGICAL             FIRST

C## S A V E:
                             SAVE FIRST, LASTPR, POS1, POS2, LSTNAM

C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A
                             DATA FIRST/T/, LASTPR/0/, LSTNAM/BLANK/

C##                                                E X E C U T I O N
C##                                                E X E C U T I O N

      IF ( PNO .NE. LASTPR .AND. LASTPR .NE. 0 ) THEN
          WRITE ( UNIT, 99999 ) COMPTD
          WRITE ( UNIT, 99999 ) OUT(1)(1:POS1)
          WRITE ( UNIT, 99999 ) OUT(2)(1:POS2)
          LASTPR = 0
      ENDIF
      IF ( NAME .NE. LSTNAM ) THEN
          IF (.NOT. FIRST ) THEN
            ENDFILE UNIT
            CLOSE (UNIT)
          ENDIF
          FIRST = F
          LSTNAM = NAME
          FNAME = NAME
          IF ( NAME .NE. BLANK ) THEN
             CALL ZZCASE(FNAME,CTOUPP)
             POS = ZZLENG(FNAME) + 1
             FNAME(POS:) = TEX
             OPEN ( UNIT, FILE=PATH//FNAME, STATUS='UNKNOWN',
     -         ACCESS='SEQUENTIAL')
             REWIND UNIT
             WRITE(UNIT,*) BLANK
          ENDIF
      ENDIF

      IF ( ERR .EQ. 0 ) THEN

          OUT(1) = OBRACE
          POS = 2
          POS = POS + ZZLFTI(OUT(1)(POS:),PNO)
          OUT(1)(POS:POS+1) = CBRACE//OBRACE
          POS = POS + 2

          POS = POS + ZZLFTI(OUT(1)(POS:),ITERS)
          OUT(1)(POS:POS+1) = CBRACE//OBRACE
          POS = POS + 2

          POS = POS + ZZLFTI(OUT(1)(POS:),FVALS)
          OUT(1)(POS:POS+1) = CBRACE//OBRACE
          POS = POS + 2

          POS = POS + ZZLFTI(OUT(1)(POS:),ACC)
          IF ( CRIT .NE. 1 ) THEN
              OUT(1)(POS:POS+1) = COMMA//BLANK
              POS = POS + 2
              POS = POS + ZZLFTI(OUT(1)(POS:),CRIT)
          ENDIF
          OUT(1)(POS:POS+1) = CBRACE//OBRACE
          WRITE(OUT(1)(POS+2:),99998)
     -      MTIME,BSLASH,BSLASH,BSLASH,BSLASH,
     -      MAX(1,MIN(99,NINT(FTIME*100.0/MTIME))),BSLASH
C     -      MTIME,BSLASH,BSLASH,BSLASH,FTIME*100.0/MTIME,BSLASH
          POS = ZZLENG(OUT(1)) + 1
          OUT(1)(POS:POS) = CBRACE
          POS1 = POS

          POS = 1
          OUT(2)(POS:) = OBRACE//FP
          POS = POS + 4
          IF ( FXSTAR .NE. ZERO ) THEN
              WRITE ( OUT(2)(POS:), '(E23.15)') FC
          ELSE
              WRITE ( OUT(2)(POS:), '(E10.2)') FC
          ENDIF
          POS = ZZLENG(OUT(2)) + 1
          OUT(2)(POS:POS+2) = BLANK//CBRACE//OBRACE
          POS = POS + 3

          OUT(2)(POS:) = FP
          POS = POS + 3
          WRITE ( OUT(2)(POS:), '(E10.2)') GNORM
          POS = ZZLENG(OUT(2)) + 1
          OUT(2)(POS:POS+2) = BLANK//CBRACE//OBRACE
          POS = POS + 3

          OUT(2)(POS:) = FP
          POS = POS + 3
          WRITE ( OUT(2)(POS:), '(E10.2)') GMIN
          POS = ZZLENG(OUT(2)) + 1
          OUT(2)(POS:POS+2) = BLANK//CBRACE//OBRACE
          POS = POS + 3

          IF ( DISTX .GE. ZERO ) THEN
              OUT(2)(POS:) = FP
              POS = POS + 3
              WRITE ( OUT(2)(POS:), '(E10.2)') DISTX
          ELSE
              OUT(2)(POS:) = NOX
          ENDIF
          POS = ZZLENG(OUT(2)) + 1
          OUT(2)(POS:POS+1) = BLANK//CBRACE
          POS2 = POS + 1

          IF ( GNORM .NE. ZERO ) THEN
             WRITE (99,*) NAME,' #',PNO,' OK     AT ACC=',ACC
             LASTPR = PNO
          ELSE
             WRITE (99,*) NAME,' #',PNO,' OK     AT ACC=',ACC,
     -          ' GNORM= ',GNORM
             WRITE (99,*) ' '
          ENDIF
      ELSE IF ( NAME .NE. BLANK ) THEN
          WRITE (99,*) NAME,' #',PNO,' FAILED AT ACC=',ACC, ',ERR= ',ERR
          WRITE (99,*) ' '
      ENDIF

C## E X I T
90000      RETURN

C## F O R M A T S:

C99998 FORMAT (F7.2,A1,'SMALL',A1,'SL',A1,'S1(',F5.1,A1,'%)')
99998 FORMAT (F8.2,A1,'SMALL',A1,'SL',A1,'S1(',A1,'FIXTWO{'
     -                            ,I2,'}',A1,'%)')
99999 FORMAT ( A )

C##                 E N D OF ZZTEX.
                    END
