
      INTEGER FUNCTION ITOA(NUMBER,STRING)
*
* converts an integer to an equivalent string. The length of the
* string is returned through the function name
*     .. Parameters ..
      CHARACTER SPACE
      PARAMETER (SPACE=' ')
*     ..
*     .. Scalar Arguments ..
      INTEGER NUMBER
      CHARACTER STRING* (*)
*     ..
*     .. Local Scalars ..
      INTEGER ALEN,DIGIT,I,J,SLEN,TEMP
      CHARACTER CH,NUMBS* (10)
*     ..
*     .. External Subroutines ..
      EXTERNAL ERROR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC LEN,MOD
*     ..
*     .. Data statements ..
      DATA NUMBS/'0123456789'/
*     ..
*     .. Executable Statements ..
* initialize
      ALEN = 0
      SLEN = LEN(STRING)
      DO 10 I = 1,SLEN
          STRING(I:I) = SPACE
   10 CONTINUE
* repeat loop until temp=0
      TEMP = NUMBER
   20 DIGIT = MOD(TEMP,10) + 1
      TEMP = TEMP/10
      ALEN = ALEN + 1
      IF (ALEN.GT.SLEN) THEN
          CALL ERROR('String not long enough - in itoa')
      END IF

      STRING(ALEN:ALEN) = NUMBS(DIGIT:DIGIT)
      IF (TEMP.NE.0) GO TO 20
* reverse the string
      J = ALEN
      DO 30 I = 1,ALEN/2
          CH = STRING(I:I)
          STRING(I:I) = STRING(J:J)
          STRING(J:J) = CH
          J = J - 1
   30 CONTINUE
      ITOA = ALEN
      END
      INTEGER FUNCTION LENARR(ARR,MAXARR)
*
*.. find the 'length' of a character array i.e., the position
*.. of the first non-space character from the end.
*     .. Parameters ..
      CHARACTER SPACE
      PARAMETER (SPACE=' ')
*     ..
*     .. Scalar Arguments ..
      INTEGER MAXARR
*     ..
*     .. Array Arguments ..
      CHARACTER ARR(MAXARR)
*     ..
*     .. Local Scalars ..
      INTEGER I
*     ..
*     .. Executable Statements ..
      DO 10 I = MAXARR,1,-1
          IF (ARR(I).NE.SPACE) GO TO 20
   10 CONTINUE
      I = 0
   20 LENARR = I
      END
      SUBROUTINE OUTLIN(LINE,LENGTH)
*
* routine to output line(1:length) to the output file
* the associated unit number is nout in the common block
* cbib (see outset)
*     .. Scalar Arguments ..
      INTEGER LENGTH
      CHARACTER LINE* (*)
*     ..
*     .. Scalars in Common ..
      INTEGER NOUT
*     ..
*     .. Common blocks ..
      COMMON /CBIB/NOUT
*     ..
*     .. Save statement ..
      SAVE /CBIB/
*     ..
*     .. Executable Statements ..
      WRITE (NOUT,FMT='(a)') LINE(1:LENGTH)
      END
      SUBROUTINE SETUP
*
* ***** Implementation Dependent Routine *****
*
* nin - unit number for data reference file
* nerr - unit number for error messages
* filenm - implementation dependent name of data reference file
*
*     .. Scalars in Common ..
      INTEGER NERR,NIN,TYPE
      LOGICAL EOF
*     ..
*     .. Local Scalars ..
      CHARACTER FILENM*15
*     ..
*     .. External Subroutines ..
      EXTERNAL ERROR
*     ..
*     .. Common blocks ..
      COMMON /CIO/EOF,TYPE,NIN,NERR
*     ..
*     .. Save statement ..
      SAVE /CIO/
*     ..
*     .. Executable Statements ..
      NIN = 5
      NERR = 6
      FILENM = 'acm.dat'
*
* ***** Open statement may be implementation dependent *****
      OPEN (NIN,FILE=FILENM,STATUS='old',ACCESS='sequential',
     +     FORM='formatted',ERR=10)
      RETURN
*
* deal with error condition on opening file
   10 CALL ERROR('Failed to open data file')
      END
