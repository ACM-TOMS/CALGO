      SUBROUTINE ZZSORT ( SORT, LIST, LISTNO, PRECNO, NAMES, KLEN, ASC )

C## A R G U M E N T S:
                      INTEGER        SORT, KLEN
                      INTEGER        LIST(*),  LISTNO,  PRECNO(3,*)
                      LOGICAL        ASC
                      CHARACTER*(*)  NAMES

C## S T A T U S:
C               SINGLE/DOUBLE CONVERSION:        NOT REQUIRED.
C               SYSTEM  DEPENDENCE:                      NONE.
C
C>RCS $HEADER: SORT.F,V 1.10 91/11/20 10:53:08 BUCKLEY EXP $
C>RCS $LOG:     SORT.F,V $
C>RCS REVISION 1.10  91/11/20  10:53:08  BUCKLEY
C>RCS FINAL SUBMISSION TO TOMS
C>RCS
C>RCS REVISION 1.9  89/06/30  13:39:52  BUCKLEY
C>RCS PREPARING SUBMITTED VERSION OF MT
C>RCS
C>RCS REVISION 1.3.1.1  89/05/20  16:42:59  BUCKLEY
C>RCS TEMP. TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.3  89/05/18  12:20:52  BUCKLEY
C>RCS FINAL TEST OF MT BEFORE SUBMITTING
C>RCS
C>RCS REVISION 1.2  89/05/15  14:48:16  BUCKLEY
C>RCS INITIAL INSTALLATION OF MT INTO RCS FORM.
C>RCS
C>RCS REVISION 1.1  89/01/17  16:50:13  BUCKLEY
C>RCS INITIAL REVISION
C>RCS
C
C## D E S C R I P T I O N:
C
C     THIS ROUTINE SORTS THE LIST OF PROBLEMS, NAMELY LIST, ACCORDING
C     TO ONE OF SEVERAL POSSIBILITIES.  THE METHOD IS BASICALLY THAT
C     OF QUICK SORT DUE TO HOARE BUT AS MODIFIED BY RICHARD
C     C. SINGLETON ( REF. ALGORITHM 347, CACM ) .
C
C     THERE IS THEREFORE A LIMIT TO THE MAXIMUM SIZE OF ARRAY WHICH
C     CAN BE SORTED.  SPECIFICALLY, LISTNO IS THE NUMBER OF ELEMENTS
C     IN THE LIST TO BE SORTED, AND ON ENTRY YOU MUST HAVE
C
C               LISTNO < 2**(SIZE+1).
C
C     IN THE CURRENT VERSION, SIZE = 12, SO UP TO 2**13 - 1= 8095
C     ENTRIES MAY BE IN THE LIST TO SORT. THAT SHOULD BE ENOUGH
C     FOR HERE! SIZE IS A PARAMETER GIVEN BELOW.
C
C     THE ARGUMENTS ARE THE FOLLOWING:
C
C        SORT    SORT CODE; SEE BELOW.
C
C        LIST    THE ARRAY OF ENTRIES TO BE SORTED.
C        LISTNO  THE NUMBER OF ELEMENTS IN LIST.
C
C        PRECNO  AN ARRAY GIVING THE RECORD NUMBER OF EACH PROBLEM IN
C                THE DIRECT ACCESS FILE. THERE ARE THREE ENTRIES FOR
C                EACH PROBLEM, ITS RECORD NUMBER, ITS DIMENSION AND
C                THE FUNCTION NUMBER USED BY THE PROBLEM.
C        NAMES   A CHARACTER STRING GIVING THE NAME OF EACH PROBLEM
C                AND THE NAME OF THE FUNCTION USED BY THE PROBLEM.
C        KLEN    THE LENGTH OF EACH NAME IN THE STRING NAMES.
C        ASC     THE SORT IS IN ASCENDING ORDER IF THIS FLAG IS TRUE;
C                OTHERWISE IT IS DESCENDING.
C
C     THE SORT CODES ARE THE FOLLOWING.
C
C        SPRNAM  BY PROBLEM NAME.
C        SPRNUM  BY PROBLEM NUMBER.
C        SRECNO  BY RECORD NUMBER OF THE PROBLEM IN THE DAUF FILE
C                   (SAME AS ORDER IN THE PROLOG FILE).
C        SDIMN   BY THE (FIRST) DIMENSION OF THE PROBLEM.
C        SFNNAM  BY THE NAME OF THE FUNCTION USED WITH EACH PROBLEM.
C        SFNNUM  BY THE NUMBER OF THE FUNCTION USED WITH EACH PROBLEM.
C        SASIS   NO SORT; THE ORDER IS LEFT AS IS.
C
C## E N T R Y   P O I N T S: THE NATURAL ENTRY ZZSORT.
C
C## S U B R O U T I N E S:
C
C            LLT ...   INTRINSIC
C            LOC       A STATEMENT FUNCTION LOCATION IN NAMES ARRAY.
C            COMPRA    A STATEMENT FUNCTION FOR COMPARISON.
C            COMPRD    A STATEMENT FUNCTION FOR COMPARISON.
C            COMPAR    A STATEMENT FUNCTION FOR COMPARISON.
C
C## P A R A M E T E R S:

C-----                      TO LIMIT THE SIZE OF SORTS.
      INTEGER     SIZE
      PARAMETER ( SIZE = 12 )

C-----                      DEFINE THE POSITIONS IN PRECNO.
C   DEFINITIONS OF THE ROWS IN THE ARRAY PRECNO. ( PRECNO HOLDS THE
C   RECORD NUMBER IN THE DAUF FILE, THE MINIMUM DIMENSION, AND THE
C   FUNCTION NUMBER OF EACH PROBLEM. )

      INTEGER     RECN,     DIMN,     FNO1
      PARAMETER ( RECN = 1, DIMN = 2, FNO1 = 3 )

C-----                      DEFINE THE TYPES OF SORTS.

      INTEGER     INDEX,     NUMERC,     CHARAC
      PARAMETER ( INDEX = 1, NUMERC = 2, CHARAC = 3 )

C-----                      DEFINE CODES FOR SORTING

      INTEGER     SPRNAM,     SPRNUM,     SASIS,     SRECNO
      PARAMETER ( SPRNAM = 1, SPRNUM = 2, SASIS = 3, SRECNO = 4 )

      INTEGER     SFNNAM,     SFNNUM,     SDIMN,     SPROLG
      PARAMETER ( SFNNAM = 5, SFNNUM = 6, SDIMN = 7, SPROLG = 8 )

C## L O C A L   D E C L:

C-----LOCAL VARIABLES.

      INTEGER   I, II, IJ,  J,  K,  L,  M,  PROBI, PROBJ,  PROBK, ROW
      INTEGER   PROBL,  PROBT,  PROBTT,  TEMP, TYPE, K1, LOC, BASE

      LOGICAL        COMPAR, COMPRA, COMPRD

C-----THE DIMENSION OF THE FOLLOWING MAY NEED UPPER ADJUSTMENT IF
C-----THE NUMBER OF PROBLEMS INCREASES. SEE THE ABOVE DESCRIPTION
C-----AND THE PARAMETER SIZE.

      INTEGER        IL(0:SIZE-1),  IU(0:SIZE-1)

C## S A V E:
            SAVE  ROW, BASE, TYPE, K1

C## E Q U I V A L E N C E S: NONE ARE DEFINED.
C## C O M M O N:             NONE IS DEFINED.
C## D A T A:                 NONE ARE SET.
C##                                                E X E C U T I O N
C##                                                E X E C U T I O N
C---------          STATEMENT FUNCTION DEFINITIONS.

      LOC(I) = BASE + (I-1)*2*KLEN + 1

      COMPRA  (I,J) =
     -          TYPE .EQ. NUMERC .AND. PRECNO(ROW,I) .LT. PRECNO(ROW,J)
     -     .OR. TYPE .EQ. INDEX  .AND.     I         .LT.    J
     -     .OR. TYPE .EQ. CHARAC .AND.  LLT( NAMES(LOC(I):LOC(I)+K1),
     -                                       NAMES(LOC(J):LOC(J)+K1) )

      COMPRD  (I,J) =
     -          TYPE .EQ. NUMERC .AND. PRECNO(ROW,I) .GT. PRECNO(ROW,J)
     -     .OR. TYPE .EQ. INDEX  .AND.     I         .GT.    J
     -     .OR. TYPE .EQ. CHARAC .AND.  LGT( NAMES(LOC(I):LOC(I)+K1),
     -                                       NAMES(LOC(J):LOC(J)+K1) )

      COMPAR(I,J)=ASC .AND. COMPRA(I,J) .OR. .NOT. ASC .AND. COMPRD(I,J)

C-----IF THE SORT IS AS IS, DO NOTHING.

      IF ( SORT .EQ. SASIS ) THEN
         GO TO 90000
      ELSE
         IF ( LISTNO .GE. 2**(SIZE+1) ) THEN
            STOP 'TOO MANY TO SORT! IN ZZSORT'
         ENDIF

         K1   = KLEN - 1
         BASE = 0
         ROW  = 1

         IF      ( SORT .EQ. SPRNAM ) THEN
            BASE = 0
            TYPE = CHARAC
         ELSE IF ( SORT .EQ. SFNNAM ) THEN
            BASE = KLEN
            TYPE = CHARAC
         ELSE IF ( SORT .EQ. SPRNUM ) THEN
            TYPE = INDEX
         ELSE IF ( SORT .EQ. SRECNO .OR. SORT .EQ. SPROLG ) THEN
            TYPE = NUMERC
            ROW  = RECN
         ELSE IF ( SORT .EQ. SFNNUM ) THEN
            TYPE = NUMERC
            ROW  = FNO1
         ELSE IF ( SORT .EQ. SDIMN ) THEN
            TYPE = NUMERC
            ROW  = DIMN
         ENDIF
      ENDIF

C-----OTHERWISE, SORT IN THE MANNER DESIRED.

      M  = 0
      I  = 1
      II = 1
      J  = LISTNO
      GO TO 400

  100 IJ = ( I + J ) / 2
      K  =   I
      L  =       J

      PROBT = LIST(IJ)
      PROBI = LIST(I )
      PROBJ = LIST( J)

      IF ( COMPAR ( PROBT, PROBI ) ) THEN
         LIST(IJ) = PROBI
         LIST(I ) = PROBT
         TEMP     = PROBT
         PROBT    = PROBI
         PROBI    = TEMP
      ENDIF

      IF ( COMPAR ( PROBJ, PROBT ) ) THEN
         LIST(IJ) = PROBJ
         LIST( J) = PROBT
         TEMP     = PROBT
         PROBT    = PROBJ
         PROBJ    = TEMP

         IF ( COMPAR ( PROBT, PROBI ) ) THEN
            LIST(IJ) = PROBI
            LIST(I ) = PROBT
            TEMP     = PROBT
            PROBT    = PROBI
            PROBI    = TEMP
         ENDIF

      ENDIF

  200 L     = L - 1
      PROBL = LIST(L)

      IF ( COMPAR ( PROBT, PROBL ) ) THEN
         GO TO 200
      ELSE
         PROBTT = PROBL
      ENDIF

  300 K     = K + 1
      PROBK = LIST(K)

      IF     (  COMPAR ( PROBK, PROBT )  ) THEN
         GO TO 300
      ELSEIF ( K .LE. L ) THEN
         LIST(L) = LIST(K)
         LIST(K) = PROBTT
         GO TO 200
      ELSEIF ( L - I .GT. J - K ) THEN
         IL(M) = I
         IU(M) = L
         I     = K
      ELSE
         IL(M) = K
         IU(M) = J
         J     = L
      ENDIF

      M = M + 1

  400 IF     ( J - I .GT. 10 ) THEN
         GO TO 100
      ELSEIF (     I .EQ. II ) THEN
         IF ( I .LT. J ) THEN
            GO TO 100
         ENDIF
      ENDIF
            TEMP = I + 1
            DO 420 I = TEMP, J
               PROBT = LIST(I)
               K     = I - 1
               PROBK = LIST(K)

               IF ( COMPAR ( PROBT, PROBK ) ) THEN
  410             LIST(K + 1) = LIST(K)
                  K           = K - 1
                  PROBK       = LIST(K)
                  IF ( COMPAR ( PROBT, PROBK ) ) THEN
                     GO TO 410
                  ELSE
                     LIST(K + 1) = PROBT
                  ENDIF
               ENDIF
  420       CONTINUE

            M = M - 1
            IF ( M .GE. 0 ) THEN
               I = IL(M)
               J = IU(M)
               GO TO 400
            ENDIF
      GO TO 90000

C## E X I T
90000       RETURN

C## F O R M A T S:  NONE ARE DEFINED.
C##                 E N D         OF ZZSORT.
                    END
