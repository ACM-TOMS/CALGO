
      PROGRAM BIBTEX
*
* program to generate a bibtex data base from the data
* reference file
*
*     .. Parameters ..
      INTEGER MAXBUF
      PARAMETER (MAXBUF=80)
      INTEGER AUTLEN,TITLEN,KEYLEN,OTHLEN
      PARAMETER (AUTLEN=80,TITLEN=160,KEYLEN=400,OTHLEN=300)
      INTEGER NCLASS
      PARAMETER (NCLASS=57)
      INTEGER MAXALG
      PARAMETER (MAXALG=750)
*     ..
*     .. Scalars in Common ..
      INTEGER MAXSHR,NERR,NIN,NOALGS,NOUT,NUMBER,PAGEND,PAGEST,TYPE,
     +        VOLUME,YEAR
      LOGICAL EOF
      CHARACTER BUFFER* (MAXBUF),ALABEL* (6),JOURNL* (4),MONTH* (9),
     +          LANG* (3),SHARE* (3)
*     ..
*     .. Arrays in Common ..
      INTEGER LENREF(MAXALG),LENTIT(MAXALG),PAGE(MAXALG),
     +        SHARPT(NCLASS+1),SHSORT(MAXALG),VOL(MAXALG)
      CHARACTER AUTHOR(AUTLEN),JOUR(MAXALG),KEYWDS(KEYLEN),
     +          OTHERS(OTHLEN),REFS(OTHLEN,MAXALG),TITLE(TITLEN),
     +          TITLES(TITLEN,MAXALG),SHAR(MAXALG)*3,LABEL(MAXALG)*6,
     +          SHARMN(NCLASS)* (3),SHARST(NCLASS)* (45)
*     ..
*     .. Common blocks ..
*
* common block for storing information required to
* build an accumulative index based on the shar
* classification
*
* For each of the noalgs algorithms
* titles (character array) stores full titles,
* refs (character array) stores the comments and remarks
*                          in compressed form
* jour   (character array) holds the journal code
*                          (c=cacm, t=toms, x=others)
* shar   (character array) holds the shar index
* vol    (integer array)   holds the volume number
* page   (integer array)   holds the start page of the article
* label  (character array) stores the algorithm number
*
* lentit (integer array)   stores the number of significant
*                          characters in the title
* lenoth (integer array)   stores the number of significant
*                          characters in the other entries
* shsort (integer array)   contains a sorted list of indexes
*                          by shar classification
*
      COMMON /CALLIN/VOL,PAGE,NOALGS,LENTIT,LENREF,SHSORT
      COMMON /CALLST/TITLES,REFS,JOUR,SHAR,LABEL
*
* nout is used to define the output channel for
* the bibtex database -- see subroutine outcum
*

      COMMON /CBIB/NOUT
*
*  eof    (logical) - true if no more references to process
*                     false otherwise
*  type   (integer) - type of current input line
*                     newfld (1) - a new field
*                     cont   (2) - a continuation line
*                     newalg (3) - a new reference entry
*                     newfld, cont and newalg are parameter values
*  nin    (integer) - unit number for data reference file
*  nerr   (integer) - unit number for error messages
* 
      COMMON /CIO/EOF,TYPE,NIN,NERR
*
*  buffer (character*(maxbuf)) - used to buffer each line of data
*                                prior to processing. For the 
*                                distributed version of the data
*                                file maxbuf may be set to 80.
*
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
*
* common block containing information on shar classification index
*
* sharmn (character array) containing the shar mnemonics
* sharst (character array) containing explanatory strings
* sharpt (integer array)   used to store pointers into the
*                          array containing sorted list of
*                          algorithms
*
* maxshr (integer) maximum number of shar classifications
*
      COMMON /CSHAR/SHARMN,SHARST
      COMMON /CSHARN/MAXSHR,SHARPT
*     ..
*     .. Save statement ..
      SAVE /CBIB/,/CIO/,/CIOSTR/,/CREFNO/,/CREFST/,/CSHAR/,/CSHARN/,
     +     /CALLST/,/CALLIN/

*     .. External Functions ..
      LOGICAL GETREF
      EXTERNAL GETREF
*     ..
*     .. External Subroutines ..
      EXTERNAL INITRF,OUTREF,OUTSET,PRMBLE,SETUP
*     ..
*     .. Executable Statements ..
* set up machine dependent values
      CALL SETUP
* initialize data file system
      CALL INITRF
* open output file
      CALL OUTSET
      CALL PRMBLE
* keep processing until the end of file
   10 IF (GETREF()) THEN
          CALL OUTREF
          GO TO 10
      END IF
      END
      INTEGER FUNCTION INDARR(ARRAY,LENARR,ABEG,SPTSTR)
*
*.. return the element of the start of the string sptstr in
*.. character array array; zero if not present
*
*     .. Scalar Arguments ..
      INTEGER ABEG,LENARR
      CHARACTER SPTSTR* (*)
*     ..
*     .. Array Arguments ..
      CHARACTER ARRAY(LENARR)
*     ..
*     .. Local Scalars ..
      INTEGER I,J,SPTLEN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC LEN
*     ..
*     .. Executable Statements ..
      SPTLEN = LEN(SPTSTR)
      DO 20 I = ABEG - 1,LENARR - SPTLEN
          DO 10 J = 1,SPTLEN
              IF (ARRAY(I+J).NE.SPTSTR(J:J)) GO TO 20
   10     CONTINUE
*.. match
          INDARR = I + 1
          RETURN
   20 CONTINUE
*.. no match
      INDARR = 0
      END
      SUBROUTINE OUTAUT(AUTHOR,AUTLEN)
*
* output the author field with quotes in the form
*     author = string
* author is a string containing possibly multiple
* authors separated by semi colons.
* The output field will be restricted to lines of
* width maxlin split after either an authors name
* or an " and ' separator.
* land is the length of the string andstr
*
*     .. Parameters ..
      INTEGER MAXLIN,MAXAUT,LAND
      PARAMETER (MAXLIN=80,MAXAUT=10,LAND=5)
      CHARACTER*(*) QUOTE,COMMA,SEMIC,ANDSTR
      PARAMETER (QUOTE='"',COMMA=',',SEMIC=';',ANDSTR=' and ')
*     ..
*     .. Scalar Arguments ..
      INTEGER AUTLEN
*     ..
*     .. Array Arguments ..
      CHARACTER AUTHOR(AUTLEN)
*     ..
*     .. Local Scalars ..
      INTEGER ALEN,IARR,IST,NARRAY
      CHARACTER LINE* (MAXLIN)
*     ..
*     .. Local Arrays ..
      CHARACTER AUTARR(MAXAUT)* (MAXLIN)
*     ..
*     .. External Functions ..
      INTEGER LENARR,LENGTH
      EXTERNAL LENARR,LENGTH
*     ..
*     .. External Subroutines ..
      EXTERNAL ERROR,OUTLIN,SPLIT
*     ..
*     .. Executable Statements ..
*
* split the line into its narray authors
      ALEN = LENARR(AUTHOR,AUTLEN)
      CALL SPLIT(AUTHOR,ALEN,AUTARR,NARRAY,SEMIC)
* initialize
      LINE(1:12) = '    author='//QUOTE
      IST = 12
      IARR = 1
* loop through the authors
   10 ALEN = LENGTH(AUTARR(IARR))
* don't split an authors name
      IF (ALEN.GT.MAXLIN-1) THEN
          CALL ERROR('Increase maxlin in outaut')
      END IF
* see if there is room for this author
      IF (IST+ALEN.LT.MAXLIN) THEN
          LINE(IST+1:IST+ALEN) = AUTARR(IARR) (1:ALEN)
          IST = IST + ALEN
      ELSE
* if not flush current line
          CALL OUTLIN(LINE,IST)
          LINE(1:ALEN) = AUTARR(IARR) (1:ALEN)
          IST = ALEN
      END IF
      IF (IARR.NE.NARRAY) THEN
* more authors to come -- put the and in
          IF (IST+LAND.LT.MAXLIN) THEN
              LINE(IST+1:IST+LAND) = ANDSTR
              IST = IST + LAND
          ELSE
* not enough room
              CALL OUTLIN(LINE,IST)
              LINE(1:LAND) = ANDSTR
              IST = LAND
          END IF
          IARR = IARR + 1
          GO TO 10
      END IF
* finished - there is always room for the final quote and comma
      LINE(IST+1:IST+2) = QUOTE//COMMA
      CALL OUTLIN(LINE,IST+2)
      END
      SUBROUTINE OUTCLS
*
* output the final closing brace
*
*     .. External Subroutines ..
      EXTERNAL OUTLIN
*     ..
*     .. Executable Statements ..
      CALL OUTLIN('}',1)
      END
      SUBROUTINE OUTJRL(JOURNL)
*
* outputs the final field in the reference in the form
*      journal = cacm
* or   journal = toms
* here we are using string substitution (see prmble for
* the string definitions).
*
* Because it is the last field there is no trailing comma.
*
*     .. Parameters ..
      INTEGER MAXLIN
      PARAMETER (MAXLIN=80)
*     ..
*     .. Scalar Arguments ..
      CHARACTER JOURNL* (*)
*     ..
*     .. Local Scalars ..
      CHARACTER LINE* (MAXLIN)
*     ..
*     .. External Subroutines ..
      EXTERNAL OUTLIN
*     ..
*     .. Executable Statements ..
      LINE(1:16) = '    journal='//JOURNL
      CALL OUTLIN(LINE,16)
      END
      SUBROUTINE OUTKEY(ALABEL)
*
* output the initial line for a bibtex reference
*          @article{keyword
* The keyword is made up of the parameter string prefix and
* the label given to the algorithm on publication -- in
* all cases bar two this was an integer, in the two odd cases
* an integer with an 'a' appended.
*
*     .. Parameters ..
      INTEGER MAXLIN
      CHARACTER*(*) PREFIX,COMMA
      PARAMETER (MAXLIN=80,PREFIX='acmalg',COMMA=',')
*     ..
*     .. Scalar Arguments ..
      CHARACTER ALABEL* (*)
*     ..
*     .. Local Scalars ..
      INTEGER ALEN,LENG
      CHARACTER LINE* (MAXLIN)
*     ..
*     .. External Functions ..
      INTEGER LENGTH
      EXTERNAL LENGTH
*     ..
*     .. External Subroutines ..
      EXTERNAL ERROR,OUTLIN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC LEN
*     ..
*     .. Executable Statements ..
*
* want to strip trailing blanks
      ALEN = LENGTH(ALABEL)
* error if label field is blank
      IF (ALEN.EQ.0) THEN
          CALL ERROR('Blank label field in outkey')
      END IF
      LINE = '@article{'//PREFIX//ALABEL(1:ALEN)//COMMA
      LENG = 9 + LEN(PREFIX) + ALEN + 1
      CALL OUTLIN(LINE,LENG)
      END
      SUBROUTINE OUTMON(MONTH)
*
* output the month field in the format
*     month = string
* This routine outputs an explicit string. Lamport
* suggests using string definitions for this. The required
* changes are straightforward but the string definitions
* need to be included within the resulting .bib file
*(see comments in prmble)
*
*     .. Parameters ..
      INTEGER MAXLIN
      CHARACTER QUOTE,SPACE,COMMA
      PARAMETER (MAXLIN=80,QUOTE='"',SPACE=' ',COMMA=',')
*     ..
*     .. Scalar Arguments ..
      CHARACTER MONTH* (*)
*     ..
*     .. Local Scalars ..
      INTEGER ALEN,I,IST
      CHARACTER LINE* (MAXLIN)
*     ..
*     .. External Functions ..
      INTEGER LENGTH
      EXTERNAL LENGTH
*     ..
*     .. External Subroutines ..
      EXTERNAL OUTLIN
*     ..
*     .. Executable Statements ..
*
* initialize
      LINE(1:11) = '    month='//QUOTE
      IST = 11
*strip leading blanks from month
      ALEN = LENGTH(MONTH)
      DO 10 I = 1,ALEN
          IF (MONTH(I:I).NE.SPACE) GO TO 15
   10 CONTINUE
   15 LINE(IST+1:IST-I+ALEN+3) = MONTH(I:ALEN)//QUOTE//COMMA
      CALL OUTLIN(LINE,IST-I+ALEN+3)
      END
      SUBROUTINE OUTNUM(NUMBER,FIELD)
*
* output a numerical value with a field name field
* 1.e.     field = number
* This is used to output the volume, number and year fields.
*
*     .. Parameters ..
      INTEGER MAXLIN
      CHARACTER COMMA
      INTEGER NUMLEN
      PARAMETER (MAXLIN=80,COMMA=',',NUMLEN=6)
*     ..
*     .. Scalar Arguments ..
      INTEGER NUMBER
      CHARACTER FIELD* (*)
*     ..
*     .. Local Scalars ..
      INTEGER IST,NLEN
      CHARACTER LINE* (MAXLIN),NUMSTR* (NUMLEN)
*     ..
*     .. External Functions ..
      INTEGER ITOA,LENGTH
      EXTERNAL ITOA,LENGTH
*     ..
*     .. External Subroutines ..
      EXTERNAL OUTLIN
*     ..
*     .. Executable Statements ..
*
* initialize
*
      IST = LENGTH(FIELD)
      LINE(1:IST) = FIELD(1:IST)
      NLEN = ITOA(NUMBER,NUMSTR)
      LINE(IST+1:IST+NLEN+1) = NUMSTR(1:NLEN)//COMMA
      CALL OUTLIN(LINE,IST+NLEN+1)
      END
      SUBROUTINE OUTPGS(PAGEST,PAGEND)
*
* output the pages field in the format
*     pages = "pagest--pagend"
* or  pages = "pagest"
* if there is no range
*
*     .. Parameters ..
      INTEGER MAXLIN
      CHARACTER QUOTE,COMMA
      PARAMETER (MAXLIN=80,QUOTE='"',COMMA=',')
*     ..
*     .. Scalar Arguments ..
      INTEGER PAGEND,PAGEST
*     ..
*     .. Local Scalars ..
      INTEGER IST,NLEN
      CHARACTER LINE* (MAXLIN),NUMSTR* (6)
*     ..
*     .. External Functions ..
      INTEGER ITOA
      EXTERNAL ITOA
*     ..
*     .. External Subroutines ..
      EXTERNAL OUTLIN
*     ..
*     .. Executable Statements ..
* initialize
      LINE(1:11) = '    pages='//QUOTE
      IST = 11
* turn pagest into a string
      NLEN = ITOA(PAGEST,NUMSTR)
      LINE(IST+1:IST+NLEN) = NUMSTR(1:NLEN)
      IST = IST + NLEN
* deal with possible range
* pagend=0 if just a single page reference
      IF (PAGEND.NE.0) THEN
          LINE(IST+1:IST+2) = '--'
          IST = IST + 2
          NLEN = ITOA(PAGEND,NUMSTR)
          LINE(IST+1:IST+NLEN) = NUMSTR(1:NLEN)
          IST = IST + NLEN
      END IF
*add the final comma
      LINE(IST+1:IST+2) = QUOTE//COMMA
      CALL OUTLIN(LINE,IST+2)
      END
      SUBROUTINE OUTREF
*
* output a single bibliographic entry
*
*     .. Parameters ..
      INTEGER AUTLEN,TITLEN,KEYLEN,OTHLEN
      PARAMETER (AUTLEN=80,TITLEN=160,KEYLEN=400,OTHLEN=300)
*     ..
*     .. Scalars in Common ..
      INTEGER NUMBER,PAGEND,PAGEST,VOLUME,YEAR
      CHARACTER ALABEL* (6),JOURNL* (4),MONTH* (9),LANG* (3),SHARE* (3)
*     ..
*     .. Arrays in Common ..
      CHARACTER AUTHOR(AUTLEN),KEYWDS(KEYLEN),OTHERS(OTHLEN),
     +          TITLE(TITLEN)
*     ..
*     .. External Subroutines ..
      EXTERNAL OUTAUT,OUTCLS,OUTJRL,OUTKEY,OUTMON,OUTNUM,OUTPGS,OUTTIT
*     ..
*     .. Common blocks ..
      COMMON /CREFNO/VOLUME,NUMBER,YEAR,PAGEST,PAGEND
      COMMON /CREFST/ALABEL,JOURNL,MONTH,LANG,SHARE,AUTHOR,TITLE,KEYWDS,
     +       OTHERS
*     ..
*     .. Save statement ..
      SAVE /CREFNO/,/CREFST/
*     ..
*     .. Executable Statements ..
*
* @article{keyword
      CALL OUTKEY(ALABEL)
* author field
      CALL OUTAUT(AUTHOR,AUTLEN)
* title field
      CALL OUTTIT(TITLE,TITLEN)
* volume field
      CALL OUTNUM(VOLUME,'    volume=')
* number field
      CALL OUTNUM(NUMBER,'    number=')
* pages field
      CALL OUTPGS(PAGEST,PAGEND)
* year field
      CALL OUTNUM(YEAR,'    year=')
* month field
      CALL OUTMON(MONTH)
* journal field
      CALL OUTJRL(JOURNL)
* closing brace
      CALL OUTCLS
      END
      SUBROUTINE OUTSET
*
* open file for output of bibtex references
* unit number nout must not clash with nin and nerr
* which are initialized in routine setup
*
*     .. Scalars in Common ..
      INTEGER NOUT
*     ..
*     .. External Subroutines ..
      EXTERNAL ERROR
*     ..
*     .. Common blocks ..
      COMMON /CBIB/NOUT
*     ..
*     .. Save statement ..
      SAVE /CBIB/
*     ..
*     .. Executable Statements ..
      NOUT = 7
      OPEN (UNIT=NOUT,FILE='acm.bib',STATUS='unknown',FORM='formatted',
     +     ACCESS='sequential',ERR=10)
      RETURN
* deal with possible error opening output file
   10 CALL ERROR('Failed to open file for bibtex database')
      END
      SUBROUTINE OUTTIT(TITLE,TITLEN)
*
* output the title field in quotes
*     title = string
* The output field will be restricted to lines of
* width maxlin split on a space
*
*     .. Parameters ..
      INTEGER MAXLIN
      CHARACTER QUOTE,COMMA,SPACE
      PARAMETER (MAXLIN=80,QUOTE='"',COMMA=',',SPACE=' ')
*     ..
*     .. Scalar Arguments ..
      INTEGER TITLEN
*     ..
*     .. Array Arguments ..
      CHARACTER TITLE(TITLEN)
*     ..
*     .. Local Scalars ..
      INTEGER ALEN,I,IEND,IPTR,IST,TLEN
      CHARACTER LINE* (MAXLIN)
*     ..
*     .. External Functions ..
      INTEGER LENARR
      EXTERNAL LENARR
*     ..
*     .. External Subroutines ..
      EXTERNAL ERROR,OUTLIN,STRARR
*     ..
*     .. Executable Statements ..
*
* initialize
      LINE(1:11) = '    title='//QUOTE
      IST = 11
      TLEN = LENARR(TITLE,TITLEN)
      ALEN = TLEN
      IEND = 1
* can we output the rest of the title on one line
   10 IF (IST+ALEN.GT.MAXLIN-2) THEN
* break on a space
          IPTR = MAXLIN - IST - 3
          DO 20 I = IPTR + IEND,IEND,-1
              IF (TITLE(I).EQ.SPACE) THEN
                  GO TO 30

              END IF

   20     CONTINUE
          CALL ERROR('Increase maxlin in outtit')
* cut here
   30     CALL STRARR(LINE,IST+1,TITLE,IEND,I-IEND-1)
          CALL OUTLIN(LINE,IST+I-IEND)
          IEND = I + 1
          ALEN = TLEN - IEND + 1
          IST = 0
          GO TO 10

      ELSE
          CALL STRARR(LINE,IST+1,TITLE,IEND,ALEN-1)
          IST = IST + ALEN
      END IF
* it fits output with closing quote and comma
      LINE(IST+1:IST+2) = QUOTE//COMMA
      CALL OUTLIN(LINE,IST+2)
      END
      SUBROUTINE PRMBLE
*
* output any preamble to the .bib file
* Here we output string definitions for cacm or toms.
* If string substitution were to be used for months
* the definitions could also be placed here.
*
*     .. Parameters ..
      INTEGER MAXLIN
      PARAMETER (MAXLIN=80)
*     ..
*     .. Local Scalars ..
      CHARACTER LINE* (MAXLIN)
*     ..
*     .. External Subroutines ..
      EXTERNAL OUTLIN
*     ..
*     .. Executable Statements ..
      LINE(1:39) = '@string{toms="ACM Trans. Math. Softw."}'
      CALL OUTLIN(LINE,39)
      LINE(1:27) = '@string{cacm="Commun. ACM"}'
      CALL OUTLIN(LINE,27)
      LINE(1:26) = '@string{topl="ACM TOPLAS"}'
      CALL OUTLIN(LINE,26)
      END
      SUBROUTINE SPLIT(STRING,LENSTR,ARRAY,NARRAY,SPTSTR)
*
* routine split splits the string string into narray pieces
* using the separator string sptstr. Substrings are stored
* in the array array.
*
*     .. Parameters ..
      CHARACTER SPACE
      PARAMETER (SPACE=' ')
*     ..
*     .. Scalar Arguments ..
      INTEGER LENSTR,NARRAY
      CHARACTER SPTSTR* (*)
*     ..
*     .. Array Arguments ..
      CHARACTER STRING(LENSTR),ARRAY(*)* (*)
*     ..
*     .. Local Scalars ..
      INTEGER I,IEND,IST,SPTLEN
*     ..
*     .. External Functions ..
      INTEGER INDARR
      EXTERNAL INDARR
*     ..
*     .. External Subroutines ..
      EXTERNAL STRARR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC LEN
*     ..
*     .. Executable Statements ..
      SPTLEN = LEN(SPTSTR)
      NARRAY = 0
      IST = 1
   10 IEND = INDARR(STRING,LENSTR,IST,SPTSTR)
      IF (IEND.NE.0) THEN
          NARRAY = NARRAY + 1
*
* lose leading spaces
          DO 20 I = IST,IEND - 1
              IF (STRING(I).NE.SPACE) THEN
                  GO TO 30
              END IF
   20     CONTINUE
   30     CONTINUE
          ARRAY(NARRAY) = SPACE
          CALL STRARR(ARRAY(NARRAY),1,STRING,I,IEND-1-I)
          IST = IEND + SPTLEN
          GO TO 10
      END IF
      END
      SUBROUTINE STRARR(STRING,SBEG,ARRAY,ABEG,LENG)
*
*.. copy array(abeg) .. array(abeg+leng) to string(sbeg:sbeg+leng)
*
*     .. Scalar Arguments ..
      INTEGER ABEG,LENG,SBEG
      CHARACTER STRING* (*)
*     ..
*     .. Array Arguments ..
      CHARACTER ARRAY(ABEG+LENG)
*     ..
*     .. Local Scalars ..
      INTEGER I
*     ..
*     .. Executable Statements ..
      DO 10 I = 0,LENG
          STRING(SBEG+I:SBEG+I) = ARRAY(ABEG+I)
   10 CONTINUE
      END

