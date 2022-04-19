
      LOGICAL FUNCTION GETREF()
*
* Getref accumulates a complete reference to a single algorithm.
* The function returns false if no reference has been found,
* i.e., we are at the end of the file.

*
*  TEMPLATE FOR USE OF GETREF
*
*  set up i/o channels and open data file
*  (This routine contains a possibly machine dependent
*   OPEN statement)
*     CALL SETUP
*
* set up output file -- application dependent routine
*     CALL OUTFIL
*
* initialize input buffer for references
* a call to initrf must precede calls to getref
*     CALL INITRF
*
* process all references
*  10 IF (GETREF()) THEN
*         process current reference
*         GO TO 10
*     END IF
*
*
* Portability:
*
* This code along with the two example progams have been checked
* with the Toolpack1/2 semantic analyzer and Fortran 77 Pfort.
* The only warning generated was for the use of non-standard
* characters (lower case - used for comments and some format
* statements and backslash used in format statements in the
* example programs only).
*
* The example programs have been run successfully on:
*
*    SUN 3 (SUNOS 4.0) f77 compiler
*  VAX8800 (VMS V4.7)  VMS Fortran
*  VAX8800 (VMS V4.7)  Watfor77
*
*
*  buffer (character*(maxbuf)) - used to buffer each line of data
*                                prior to processing. For the
*                                distributed version of the data
*                                file maxbuf may be set to 80.
* form holds the format for splitting up the initial line of the
* reference. This contains
*
*  algorithm label, journal, start page, end page, volume, number,
*  month, year, share classification index, source language
*
*     .. Parameters ..
      INTEGER NEWFLD,NEWALG
      PARAMETER (NEWFLD=1,NEWALG=3)
      INTEGER MAXBUF
      PARAMETER (MAXBUF=80)
      INTEGER AUTLEN,TITLEN,KEYLEN,OTHLEN
      PARAMETER (AUTLEN=80,TITLEN=160,KEYLEN=400,OTHLEN=300)
      CHARACTER*(*) FORM
      PARAMETER (FORM='(a6,1x,a4,2i5,2i3,1x,a9,i5,1x,a3,1x,a3)')
      CHARACTER BLNK
      PARAMETER (BLNK=' ')
*     ..
*     .. Scalars in Common ..
      INTEGER NERR,NIN,NUMBER,PAGEND,PAGEST,TYPE,VOLUME,YEAR
      LOGICAL EOF
      CHARACTER BUFFER* (MAXBUF),ALABEL* (6),JOURNL* (4),MONTH* (9),
     +          LANG* (3),SHARE* (3)
*     ..
*     .. Arrays in Common ..
      CHARACTER AUTHOR(AUTLEN),KEYWDS(KEYLEN),OTHERS(OTHLEN),
     +          TITLE(TITLEN)
*     ..
*     .. Local Scalars ..
      INTEGER I
*     ..
*     .. External Subroutines ..
      EXTERNAL ERROR,GETFLD,NXTLIN
*     ..
*     .. Common blocks ..
      COMMON /CIO/EOF,TYPE,NIN,NERR
      COMMON /CIOSTR/BUFFER
*
* common block holds reference items which may usefully be
* represented by integers.
*
*  volume - volume number
*  number - issue number (1-12 for cacm, 1-4 for toms)
*  year   - year of publication
*  pagest - start page of article
*  pagend - end page of article (set zero if whole article confined
*           to a single page
*
      COMMON /CREFNO/VOLUME,NUMBER,YEAR,PAGEST,PAGEND
*
*  stores reference fields in string form
*
*  alabel*(6)     - label used to identify an algorithm. In all cases
*                   but two this is an integer. In the odd two cases
*                   it is an integer followed by the letter 'a'.
*  journl*(4)     - either 'cacm' or 'toms'
*  month *(9)     - month of publication e.g., January, December etc.
*  lang  *(3)     - language in which algorithm is coded
*                   'F' = Fortran (66 or 77);  'A60' = Algol 60;
*                   'R' = Ratfor'           ;  'PL1' = PL1;
*  share *(3)     - Share classification index
*  author(autlen) - possibly multiple author names. All author fields
*                   end with a semi colon, multiple authors are
*                   separated by semi colons.
*  title (titlen) - possibly multiline title concatenated into a single
*                   string. Multiple blanks are compressed to a single
*                   blank.
*  keywds(keylen) - list of keywords separated by semi colons and
*                   terminated by a semi colon.
*  others(othlen) - list of remarks and certificates pertaining to the
*                   current reference separated by and terminated by a
*                   semi colon. Each entry has the following form
*
*  R or C,journal,start page,end page,volume,number,month,year,authors;
*
*                   The first field states whether the article was a
*                   Remark or Certificate.
*                   Multiple authors are separated by ' and '.
*                   All other field are as described in the main
*                   reference.
*
*  autlen, titlen, keylen, othlen are parameter values.
      COMMON /CREFST/ALABEL,JOURNL,MONTH,LANG,SHARE,AUTHOR,TITLE,KEYWDS,
     +       OTHERS
*     ..
*     .. Save statement ..
      SAVE /CIO/,/CIOSTR/,/CREFNO/,/CREFST/
*     ..
*     .. Executable Statements ..
      IF (EOF) THEN
          GETREF = .false.
          RETURN
      END IF
* buffer should contain the first line of a new reference
      IF (TYPE.NE.NEWALG) THEN
          CALL ERROR('Inconsistent data - in getref')
      END IF
* split the initial line
      READ (BUFFER,FMT=FORM) ALABEL,JOURNL,PAGEST,PAGEND,VOLUME,NUMBER,
     +  MONTH,YEAR,SHARE,LANG
      CALL NXTLIN
* get possibly multiline fields for author, title and keywords
* these are always present.
      CALL GETFLD(AUTHOR,AUTLEN)
      CALL GETFLD(TITLE,TITLEN)
      CALL GETFLD(KEYWDS,KEYLEN)
      GETREF = .true.
* could be end of file need to process this record
      IF (EOF) THEN
*
* if others field is not present set others to all blanks
* to ensure it is defined
          DO 10 I = 1,OTHLEN
              OTHERS(I) = BLNK
   10     CONTINUE
          RETURN
      END IF
      IF (TYPE.EQ.NEWFLD) THEN
          CALL GETFLD(OTHERS,OTHLEN)
      ELSE
*
* if others field is not present set others to all blanks
* to ensure it is defined
          DO 20 I = 1,OTHLEN
              OTHERS(I) = BLNK
   20     CONTINUE
      END IF
      END
      SUBROUTINE ERROR(ERRSTR)
*
* Output routine for error messages passed via errstr
* Errors are assumed fatal.
*  parameter values
*     .. Scalar Arguments ..
      CHARACTER ERRSTR* (*)
*     ..
*     .. Scalars in Common ..
      INTEGER NERR,NIN,TYPE
      LOGICAL EOF
*     ..
*     .. Common blocks ..
      COMMON /CIO/EOF,TYPE,NIN,NERR
*     ..
*     .. Save statement ..
      SAVE /CIO/
*     ..
*     .. Executable Statements ..
      WRITE (NERR,FMT='(a)') ERRSTR
      STOP
      END
      SUBROUTINE GETFLD(FLD,MAXFLD)
*
* Builds a complete field in the character variable fld.
* This may involve reading an initial line followed by
* zero, one, or more continuation lines. This routine must
* not be used to input an initial line (type=newalg) directly.
*  parameter values
*  buffer (character*(maxbuf)) - used to buffer each line of data
*                                prior to processing. For the
*                                distributed version of the data
*                                file maxbuf may be set to 80.
*     .. Parameters ..
      INTEGER NEWFLD,CONT
      PARAMETER (NEWFLD=1,CONT=2)
      INTEGER MAXBUF
      PARAMETER (MAXBUF=80)
      CHARACTER SPACE
      PARAMETER (SPACE=' ')
*     ..
*     .. Scalar Arguments ..
      INTEGER MAXFLD
*     ..
*     .. Array Arguments ..
      CHARACTER FLD(MAXFLD)
*     ..
*     .. Scalars in Common ..
      INTEGER NERR,NIN,TYPE
      LOGICAL EOF
      CHARACTER BUFFER* (MAXBUF)
*     ..
*     .. Local Scalars ..
      INTEGER I,ILEN,LENFLD
*     ..
*     .. External Functions ..
      INTEGER LENGTH
      EXTERNAL LENGTH
*     ..
*     .. External Subroutines ..
      EXTERNAL ERROR,NXTLIN
*     ..
*     .. Common blocks ..
      COMMON /CIO/EOF,TYPE,NIN,NERR
      COMMON /CIOSTR/BUFFER
*     ..
*     .. Save statement ..
      SAVE /CIO/,/CIOSTR/
*     ..
*     .. Executable Statements ..
* check for consistent data
      IF (EOF .OR. TYPE.NE.NEWFLD) THEN
          PRINT *,EOF,TYPE,NEWFLD
          CALL ERROR('Inconsistent data - in getfld')
      END IF
* initialize field - strip off leading space and trailing spaces
      DO 5 I = 1,MAXFLD
          FLD(I) = SPACE
    5 CONTINUE
      LENFLD = LENGTH(BUFFER(2:MAXBUF))
      IF (LENFLD.GT.MAXFLD) THEN
          CALL ERROR('Argument to getfld not long enough')
      ELSE
          DO 7 I = 2,LENFLD + 1
              FLD(I-1) = BUFFER(I:I)
    7     CONTINUE
      END IF
* concatenate continuation lines if present
   10 CALL NXTLIN
* could have an eof here when getting last keywords field
      IF (EOF) RETURN
      IF (TYPE.EQ.CONT) THEN
          ILEN = LENGTH(BUFFER(2:MAXBUF))
          IF (LENFLD+ILEN.GT.MAXFLD) THEN
              CALL ERROR('Argument to getfld not long enough')
          ELSE
              DO 15 I = 2,ILEN + 1
                  FLD(LENFLD+I-1) = BUFFER(I:I)
   15         CONTINUE
          END IF
          LENFLD = LENFLD + ILEN
          GO TO 10
      END IF
      END
      SUBROUTINE INITRF
*
* Routine used to initialize the input buffer. It must be
* called prior to the first call to getref and after the
* call to setup.
*  parameter values
*     .. Parameters ..
      INTEGER NEWALG
      PARAMETER (NEWALG=3)
*     ..
*     .. Scalars in Common ..
      INTEGER NERR,NIN,TYPE
      LOGICAL EOF
*     ..
*     .. External Subroutines ..
      EXTERNAL ERROR,NXTLIN
*     ..
*     .. Common blocks ..
      COMMON /CIO/EOF,TYPE,NIN,NERR
*     ..
*     .. Save statement ..
      SAVE /CIO/
*     ..
*     .. Executable Statements ..
      CALL NXTLIN
* check we have a newalg line
      IF (EOF .OR. TYPE.NE.NEWALG) THEN
          CALL ERROR('Inconsistent data file - in initrf')
      END IF
      END
      INTEGER FUNCTION LENGTH(STRING)
*
* Strips trailing spaces and returns the number of characters
* to the first non-blank from the right-hand end
*     .. Parameters ..
      CHARACTER SPACE
      PARAMETER (SPACE=' ')
*     ..
*     .. Scalar Arguments ..
      CHARACTER STRING* (*)
*     ..
*     .. Local Scalars ..
      INTEGER I,MAXLEN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC LEN
*     ..
*     .. Executable Statements ..
      MAXLEN = LEN(STRING)
      DO 10 I = MAXLEN,1,-1
          IF (STRING(I:I).NE.SPACE) GO TO 20
   10 CONTINUE
      I = 0
   20 LENGTH = I
      END
      SUBROUTINE NXTLIN
*
* Routine to read next line of input from the data file into
* buffer. If end of file is detected eof is set true (false
* otherwise). The line is classified, using type, as a new
* algorithm (type = newalg(3)), a new field (type = newfld(1))
* or a continuation line (type=cont(2)). Type values are set in
* the parameter statement associated with common block cio.
*  parameter values
*  buffer (character*(maxbuf)) - used to buffer each line of data
*                                prior to processing. For the
*                                distributed version of the data
*                                file maxbuf may be set to 80.
*     .. Parameters ..
      INTEGER NEWFLD,CONT,NEWALG
      PARAMETER (NEWFLD=1,CONT=2,NEWALG=3)
      INTEGER MAXBUF
      PARAMETER (MAXBUF=80)
      CHARACTER SPACE,PLUS
      PARAMETER (SPACE=' ',PLUS='+')
*     ..
*     .. Scalars in Common ..
      INTEGER NERR,NIN,TYPE
      LOGICAL EOF
      CHARACTER BUFFER* (MAXBUF)
*     ..
*     .. Local Scalars ..
      CHARACTER FIRSTL
*     ..
*     .. External Functions ..
      INTEGER LENGTH
      EXTERNAL LENGTH
*     ..
*     .. Common blocks ..
      COMMON /CIO/EOF,TYPE,NIN,NERR
      COMMON /CIOSTR/BUFFER
*     ..
*     .. Save statement ..
      SAVE /CIO/,/CIOSTR/
*     ..
*     .. Executable Statements ..
* get a line of text
   10 READ (NIN,FMT='(a80)',END=100) BUFFER
* ignore blank lines
      IF (LENGTH(BUFFER).EQ.0) GO TO 10
* classify line by first character
      FIRSTL = BUFFER(1:1)
      IF (FIRSTL.EQ.SPACE) THEN
          TYPE = NEWFLD
      ELSE IF (FIRSTL.EQ.PLUS) THEN
          TYPE = CONT
      ELSE
          TYPE = NEWALG
      END IF
      EOF = .false.
      RETURN
* here if end of file detected on read
  100 EOF = .true.
      END


