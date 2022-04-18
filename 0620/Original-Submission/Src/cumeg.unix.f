
      SUBROUTINE ADDSHR(SHARVL,NOALG)
*
* Add the current algorithm into the sorted shar array.
* This is pretty crude but given the program will be run
* only infrequently we think it's good enough.
*
* sharvl -- shar index associated with current algorithm
* noalg  -- index into volume, pages etc arrays for current
*           algorithm information.
*
*
*     .. Parameters ..
      CHARACTER*(*) MESS
      PARAMETER (MESS='Unrecognized share value ')
      INTEGER MAXALG,TITLEN,OTHLEN
      PARAMETER (MAXALG=750,TITLEN=160,OTHLEN=300)
      INTEGER NCLASS
      PARAMETER (NCLASS=57)
*     ..
*     .. Scalar Arguments ..
      INTEGER NOALG
      CHARACTER SHARVL* (3)
*     ..
*     .. Scalars in Common ..
      INTEGER MAXSHR,NOALGS,SVAL
*     ..
*     .. Arrays in Common ..
      INTEGER LENREF(MAXALG),LENTIT(MAXALG),PAGE(MAXALG),
     +        SHARPT(NCLASS+1),SHSORT(MAXALG),VOL(MAXALG)
      CHARACTER JOUR(MAXALG),REFS(OTHLEN,MAXALG),TITLES(TITLEN,MAXALG),
     +          SHAR(MAXALG)*3,LABEL(MAXALG)*6,SHARMN(NCLASS)* (3),
     +          SHARST(NCLASS)* (45)
*     ..
*     .. Local Scalars ..
      INTEGER I,J,POS,SHPOS
      CHARACTER SHTEMP* (3)
*     ..
*     .. External Functions ..
      INTEGER BINCHP
      EXTERNAL BINCHP
*     ..
*     .. External Subroutines ..
      EXTERNAL ERROR,SHARFX
*     ..
*     .. Common blocks ..
      COMMON /CALLIN/VOL,PAGE,NOALGS,LENTIT,LENREF,SHSORT
      COMMON /CALLST/TITLES,REFS,JOUR,SHAR,LABEL
      COMMON /CSHAR/SHARMN,SHARST
      COMMON /CSHARN/MAXSHR,SHARPT,SVAL
*     ..
*     .. Save statement ..
      SAVE /CALLST/,/CALLIN/,/CSHAR/,/CSHARN/
*     ..
*     .. Executable Statements ..
*
* find pointer position for current shar value
* converting sharvl to upper case if necessary
      SHTEMP = SHARVL
      CALL SHARFX(SHTEMP)
      POS = BINCHP(SHTEMP,SHARMN,MAXSHR)
      IF (POS.EQ.0) THEN
* stop if not found
          CALL ERROR(MESS//SHTEMP)
      END IF
*
* locate where new entry should go -- want it at
* the end of the current list to ensure that the entries
* associated with each shar classification remain sorted.
*
      SHPOS = POS
      DO 30 I = POS + 1,MAXSHR
*
* find first non-zero pointer
          IF (SHARPT(I).NE.0) THEN
              SHPOS = SHARPT(I)
*
* update this and other pointers
              DO 10 J = I,MAXSHR
                  IF (SHARPT(J).NE.0) THEN
                      SHARPT(J) = SHARPT(J) + 1
                  END IF
   10         CONTINUE
*
* insert the new entry
              DO 20 J = NOALG - 1,SHPOS,-1
                  SHSORT(J+1) = SHSORT(J)
   20         CONTINUE
              SHSORT(SHPOS) = NOALG
              IF (SHARPT(POS).EQ.0) THEN
                  SHARPT(POS) = SHPOS
              END IF
              RETURN
          END IF
   30 CONTINUE
*
* add to end of list and update pointer if necessary
      SHSORT(NOALG) = NOALG
      IF (SHARPT(SHPOS).EQ.0) THEN
          SHARPT(SHPOS) = NOALG
      END IF
      END
      BLOCK DATA BDSHAR
*
* Initialize share mnemonic and descriptive string arrays
*
*     .. Parameters ..
      INTEGER NCLASS, NCLAS1
      PARAMETER (NCLASS=57,NCLAS1=NCLASS+1)
*     ..
*     .. Scalars in Common ..
      INTEGER MAXSHR,SVAL
*     ..
*     .. Arrays in Common ..
      INTEGER SHARPT(NCLASS+1)
      CHARACTER SHARMN(NCLASS)* (3),SHARST(NCLASS)* (45)
*     ..
*     .. Common blocks ..
      COMMON /CSHAR/SHARMN,SHARST
      COMMON /CSHARN/MAXSHR,SHARPT,SVAL
*     ..
*     .. Save statement ..
      SAVE /CSHAR/,/CSHARN/
*     ..
*     .. Data statements ..
      DATA SHARMN(1)/'A1 '/
      DATA SHARST(1)/'Real Arithmetic, Number Theory'/
      DATA SHARMN(2)/'A2 '/
      DATA SHARST(2)/'Complex Arithmetic'/
      DATA SHARMN(3)/'B1 '/
      DATA SHARST(3)/'Trig and Inverse Trig Functions'/
      DATA SHARMN(4)/'B2 '/
      DATA SHARST(4)/'Hyperbolic Functions'/
      DATA SHARMN(5)/'B3 '/
      DATA SHARST(5)/'Exponential and Logarithmic Functions'/
      DATA SHARMN(6)/'B4 '/
      DATA SHARST(6)/'Roots and Powers'/
      DATA SHARMN(7)/'C1 '/
      DATA SHARST(7)/'Operations on Polynomials and Power Series'/
      DATA SHARMN(8)/'C2 '/
      DATA SHARST(8)/'Zeros of Polynomials'/
      DATA SHARMN(9)/'C5 '/
      DATA SHARST(9)/'Zeros of one or more Nonlinear Equations'/
      DATA SHARMN(10)/'C6 '/
      DATA SHARST(10)/'Summation of Series, Convergence Acceleration'/
      DATA SHARMN(11)/'D1 '/
      DATA SHARST(11)/'Quadrature'/
      DATA SHARMN(12)/'D2 '/
      DATA SHARST(12)/'Ordinary Differential Equations'/
      DATA SHARMN(13)/'D3 '/
      DATA SHARST(13)/'Partial Differential Equations'/
      DATA SHARMN(14)/'D4 '/
      DATA SHARST(14)/'Differentiation'/
      DATA SHARMN(15)/'D5 '/
      DATA SHARST(15)/'Integral Equations'/
      DATA SHARMN(16)/'E1 '/
      DATA SHARST(16)/'Interpolation'/
      DATA SHARMN(17)/'E2 '/
      DATA SHARST(17)/'Curve and Surface Fitting'/
      DATA SHARMN(18)/'E3 '/
      DATA SHARST(18)/'Smoothing'/
      DATA SHARMN(19)/'E4 '/
      DATA SHARST(19)/'Minimizing or Maximizing a Function'/
      DATA SHARMN(20)/'F1 '/
      DATA SHARST(20)/'Matrix Operations, including Inversion'/
      DATA SHARMN(21)/'F2 '/
      DATA SHARST(21)/'Eigenvalues and Eigenvectors of a Matrix'/
      DATA SHARMN(22)/'F3 '/
      DATA SHARST(22)/'Determinants'/
      DATA SHARMN(23)/'F4 '/
      DATA SHARST(23)/'Simultaneous Linear Equations'/
      DATA SHARMN(24)/'F5 '/
      DATA SHARST(24)/'Orthogonalization'/
      DATA SHARMN(25)/'G1 '/
      DATA SHARST(25)/'Simple Calculations on Statistical Data'/
      DATA SHARMN(26)/'G2 '/
      DATA SHARST(26)/'Correlation and Regression Analysis'/
      DATA SHARMN(27)/'G5 '/
      DATA SHARST(27)/'Random Number Generators'/
      DATA SHARMN(28)/'G6 '/
      DATA SHARST(28)/'Permutations and Combinations'/
      DATA SHARMN(29)/'G7 '/
      DATA SHARST(29)/'Subset Generators'/
      DATA SHARMN(30)/'H  '/
      DATA SHARST(30)/'Operations Research, Graph Structure'/
      DATA SHARMN(31)/'I5 '/
      DATA SHARST(31)/'Input -- Composite'/
      DATA SHARMN(32)/'J6 '/
      DATA SHARST(32)/'Plotting'/
      DATA SHARMN(33)/'K2 '/
      DATA SHARST(33)/'Relocation'/
      DATA SHARMN(34)/'L2 '/
      DATA SHARST(34)/'Compiling'/
      DATA SHARMN(35)/'M1 '/
      DATA SHARST(35)/'Sorting'/
      DATA SHARMN(36)/'M2 '/
      DATA SHARST(36)/'Data Conversion and Scaling'/
      DATA SHARMN(37)/'O2 '/
      DATA SHARST(37)/'Simulation of Computing Structure'/
      DATA SHARMN(38)/'R2 '/
      DATA SHARST(38)/'Symbol Manipulation'/
      DATA SHARMN(39)/'S  '/
      DATA SHARST(39)/'Approximation of Special Functions'/
      DATA SVAL/39/
      DATA SHARMN(40)/'S03'/
      DATA SHARST(40)/' '/
      DATA SHARMN(41)/'S04'/
      DATA SHARST(41)/' '/
      DATA SHARMN(42)/'S07'/
      DATA SHARST(42)/' '/
      DATA SHARMN(43)/'S09'/
      DATA SHARST(43)/' '/
      DATA SHARMN(44)/'S13'/
      DATA SHARST(44)/' '/
      DATA SHARMN(45)/'S14'/
      DATA SHARST(45)/' '/
      DATA SHARMN(46)/'S15'/
      DATA SHARST(46)/' '/
      DATA SHARMN(47)/'S16'/
      DATA SHARST(47)/' '/
      DATA SHARMN(48)/'S17'/
      DATA SHARST(48)/' '/
      DATA SHARMN(49)/'S18'/
      DATA SHARST(49)/' '/
      DATA SHARMN(50)/'S19'/
      DATA SHARST(50)/' '/
      DATA SHARMN(51)/'S20'/
      DATA SHARST(51)/' '/
      DATA SHARMN(52)/'S21'/
      DATA SHARST(52)/' '/
      DATA SHARMN(53)/'S22'/
      DATA SHARST(53)/' '/
      DATA SHARMN(54)/'S23'/
      DATA SHARST(54)/' '/
      DATA SHARMN(55)/'XX '/
      DATA SHARST(55)/'Unclassified'/
      DATA SHARMN(56)/'Y1 '/
      DATA SHARST(56)/'Physics Applications'/
      DATA SHARMN(57)/'Z  '/
      DATA SHARST(57)/'All Others'/
      DATA MAXSHR/NCLASS/
      DATA SHARPT/NCLAS1*0/
*     ..
*     .. Executable Statements ..
      END
      INTEGER FUNCTION BINCHP(TARGET,GIVEN,LENGTH)
*
* Find target in sorted array given(1..length)
* using the binary chop algorithm
*
* returns 0 if target not in given(1..length)
* returns p if given(p)=target for p=1..length
*
*     .. Scalar Arguments ..
      INTEGER LENGTH
      CHARACTER TARGET* (3)
*     ..
*     .. Array Arguments ..
      CHARACTER GIVEN(LENGTH)* (3)
*     ..
*     .. Local Scalars ..
      INTEGER L,M,U
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC LGT,LLT
*     ..
*     .. Executable Statements ..
      L = 1
      U = LENGTH
   10 IF (L.GT.U) THEN
          BINCHP = 0
          RETURN
      END IF
* halve the interval
      M = (L+U)/2
      IF (LLT(GIVEN(M),TARGET)) THEN
          L = M + 1
          GO TO 10
      ELSE IF (LGT(GIVEN(M),TARGET)) THEN
          U = M - 1
          GO TO 10
      END IF
      BINCHP = M
      END
      PROGRAM CUMIND
*
* produces latex file to generate a cumulative index of
* acm algorithms sorted by shar classification
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
      INTEGER MAXSHR,NERR,NIN,NOALGS,NOUT,NUMBER,PAGEND,PAGEST,SVAL,
     +        TYPE,VOLUME,YEAR
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
*     .. Local Scalars ..
      INTEGER CURPT,I,IPOS,LASTPT
*     ..
*     .. External Functions ..
      INTEGER GETPT
      LOGICAL GETREF
      EXTERNAL GETPT,GETREF
*     ..
*     .. External Subroutines ..
      EXTERNAL ADDSHR,INITRF,OUTALG,OUTCUM,OUTEND,OUTHED,OUTSH,SAVALG,
     +         SETUP
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
*                   'R' = Ratfor           ;  'PL1' = PL1;
*                   'L' = Lisp              ;  'F90' = Fortran 90;
*                   'M' = Matlab            ;  
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
* sval (integer) position of S section heading in the sharmn
*                and sharst arrays
*
      COMMON /CSHAR/SHARMN,SHARST
      COMMON /CSHARN/MAXSHR,SHARPT,SVAL
*     ..
*     .. Save statement ..
      SAVE /CBIB/,/CIO/,/CIOSTR/,/CREFNO/,/CREFST/,/CSHAR/,/CSHARN/,
     +     /CALLST/,/CALLIN/
*     ..
*     .. Executable Statements ..
*
*  set up i/o channels and open data file
      CALL SETUP
*
* set up output file for latex output
      CALL OUTCUM
*
* initialize input buffer for references
      CALL INITRF
      CALL OUTHED
*
* process all references
   10 IF (GETREF()) THEN
          CALL SAVALG
          GO TO 10
      END IF
*
* sort on shar index
      DO 20 I = 1,NOALGS
          J = I
*
* don't pass active do-loop variable as arg
          CALL ADDSHR(SHAR(J),J)
   20 CONTINUE
      SHARPT(MAXSHR+1) = NOALGS + 1
*
* loop through the classifications, outputting headings
* where necessary and generating entries
      IPOS = 1
      CURPT = GETPT(IPOS,MAXSHR+1,SHARPT)
   30 LASTPT = GETPT(CURPT+1,MAXSHR+1,SHARPT)
      IF (LASTPT.NE.-1) THEN
          CALL OUTSH(CURPT)
          DO 40 I = SHARPT(CURPT),SHARPT(LASTPT) - 1
              CALL OUTALG(SHSORT(I))
   40     CONTINUE
          CURPT = LASTPT
          GO TO 30
      END IF
*
* windup complete latex file
      CALL OUTEND
      END
      SUBROUTINE GETOTC(FIELD,FLEN,OTHERS,OLEN)
*
* Routine to extract a list of the citation references
* in the form e.g., C4:273 (cacm volume 4 p273).
*
* Multiple entries are separated by sepchr -- we use a
* space for the LaTeX output. No sepchr is appended to the
* end of the string.
*
* No explicit bound checks are performed on field
*
*
* field -- character array on exit holds the list of
*          references
* flen  -- number of significant characters in field
*          on exit
* others-- character array containing citation reference
* olen  -- number of significant characters in others
*
*     .. Parameters ..
      CHARACTER SEPCHR
      PARAMETER (SEPCHR=' ')
*     ..
*     .. Scalar Arguments ..
      INTEGER FLEN,OLEN
*     ..
*     .. Array Arguments ..
      CHARACTER FIELD(*),OTHERS(*)
*     ..
*     .. Local Scalars ..
      INTEGER IFLEN,IPOS
*     ..
*     .. External Functions ..
      INTEGER NXTOTH,SEMPOS
      EXTERNAL NXTOTH,SEMPOS
*     ..
*     .. Executable Statements ..
* no references
      FLEN = 0
      IF (OLEN.EQ.0) THEN
          RETURN
      END IF
* one or more
      IPOS = 1
      FLEN = 1
   10 IPOS = NXTOTH(IPOS,FIELD(FLEN),IFLEN,OTHERS)
      IPOS = SEMPOS(IPOS,OTHERS)
      FLEN = FLEN + IFLEN - 1
      IF (IPOS.NE.OLEN) THEN
          IPOS = IPOS + 1
          FIELD(FLEN+1) = SEPCHR
          FLEN = FLEN + 2
          GO TO 10
      END IF
      END
      INTEGER FUNCTION GETPT(IPOS,MAXPOS,ARRAY)
*
* Returns the first non-zero entry in the elements
* array(ipos,...,maxpos). If ipos.gt.maxpos on entry
* or all elements are zero returns -1.
*
*
*     .. Scalar Arguments ..
      INTEGER IPOS,MAXPOS
*     ..
*     .. Array Arguments ..
      INTEGER ARRAY(MAXPOS)
*     ..
*     .. Local Scalars ..
      INTEGER I,ITEMP
*     ..
*     .. Executable Statements ..
*
      ITEMP = IPOS
      IF (ITEMP.GT.MAXPOS) THEN
          GETPT = -1
      ELSE
          DO 10 I = ITEMP,MAXPOS
              IF (ARRAY(I).NE.0) THEN
                  GETPT = I
                  GO TO 20
              END IF
   10     CONTINUE
* fallen off the end
          GETPT = -1
      END IF
   20 RETURN
      END
      INTEGER FUNCTION GOTHFD(IPOS,FIELD,FLEN,OTHERS)
*
* Get next field (separated by commas) from list of remarks
* and certificates for a particular algorithm.
*
*  ipos  -- on entry points to the first character or first
*           following a comma or first following a semicolon
*           in the others array
* field -- character variable to hold field on exit
* flen  -- number of characters written into field on exit
* others-- character array holding the full other citation
*          references for a particular algorithm.
*
* This routine does not explicitly bound check field.
*
* Function returns the current position of the pointer in
* the array others (could be a comma or a semi-colon)
*
*     .. Parameters ..
      CHARACTER COMMA,SEMIC
      PARAMETER (COMMA=',',SEMIC=';')
*     ..
*     .. Scalar Arguments ..
      INTEGER FLEN,IPOS
      CHARACTER FIELD* (*)
*     ..
*     .. Array Arguments ..
      CHARACTER OTHERS(*)
*     ..
*     .. Local Scalars ..
      INTEGER ITEMP
*     ..
*     .. Executable Statements ..
      FIELD(1:1) = OTHERS(IPOS)
      ITEMP = IPOS
      FLEN = 1
   10 ITEMP = ITEMP + 1
      IF (OTHERS(ITEMP).NE.COMMA .AND. OTHERS(ITEMP).NE.SEMIC) THEN
          FLEN = FLEN + 1
          FIELD(FLEN:FLEN) = OTHERS(ITEMP)
          GO TO 10
      END IF
      GOTHFD = ITEMP
      END
      SUBROUTINE MAINRF(BUFFER,BUFLEN,NALG)
*
* construct the reference to the algorithm whose reference
* is stored in array position nalg. This will be output
* in italic in the resulting cumulative index.
* Example format: C3:49 for cacm volume 3 page 49.
*
*     .. Parameters ..
      CHARACTER SEPCHR
      PARAMETER (SEPCHR=':')
      INTEGER MAXALG,TITLEN,OTHLEN
      PARAMETER (MAXALG=750,TITLEN=160,OTHLEN=300)
*     ..
*     .. Scalar Arguments ..
      INTEGER BUFLEN,NALG
      CHARACTER BUFFER* (*)
*     ..
*     .. Scalars in Common ..
      INTEGER NOALGS
*     ..
*     .. Arrays in Common ..
      INTEGER LENREF(MAXALG),LENTIT(MAXALG),PAGE(MAXALG),SHSORT(MAXALG),
     +        VOL(MAXALG)
      CHARACTER JOUR(MAXALG),REFS(OTHLEN,MAXALG),TITLES(TITLEN,MAXALG),
     +          SHAR(MAXALG)*3,LABEL(MAXALG)*6
*     ..
*     .. Local Scalars ..
      INTEGER LENGTH
*     ..
*     .. External Functions ..
      INTEGER ITOA
      EXTERNAL ITOA
*     ..
*     .. Common blocks ..
      COMMON /CALLIN/VOL,PAGE,NOALGS,LENTIT,LENREF,SHSORT
      COMMON /CALLST/TITLES,REFS,JOUR,SHAR,LABEL
*     ..
*     .. Save statement ..
      SAVE /CALLST/,/CALLIN/
*     ..
*     .. Executable Statements ..
      BUFFER(1:1) = JOUR(NALG)
      LENGTH = ITOA(VOL(NALG),BUFFER(2:))
      BUFLEN = LENGTH + 2
      BUFFER(BUFLEN:BUFLEN) = SEPCHR
      BUFLEN = BUFLEN + 1
      LENGTH = ITOA(PAGE(NALG),BUFFER(BUFLEN:))
      BUFLEN = BUFLEN + LENGTH - 1
      END
      INTEGER FUNCTION NXTOTH(IPOS,FIELD,FLEN,OTHERS)
*
* Generate next remark or certificate reference. This involves
* extracting fields 2 (journal), 4 (volume), 3 (pages) and
* forming a string of the form C17:128. (The colon may be changed
* by altering sepchr in the parameter statement.)
*
* ipos  -- on entry points to the first character or first
*          following a semicolon in the others array
* field -- character array to hold reference on exit
* flen  -- number of characters written into field on exit
* others-- character array holding the full other citation
*          references for a particular algorithm.
*
* This routine does not explicitly bound check field.
*
* Function returns the current position of the pointer in
* the array others (could be a comma or a semi-colon)
*
*     .. Parameters ..
      CHARACTER SEPCHR,PAGSEP
      PARAMETER (SEPCHR=':',PAGSEP='-')
*     ..
*     .. Scalar Arguments ..
      INTEGER FLEN,IPOS
*     ..
*     .. Array Arguments ..
      CHARACTER FIELD(*),OTHERS(*)
*     ..
*     .. Local Scalars ..
      INTEGER I,JLEN,JPOS,PLEN,VLEN
      CHARACTER JUNK,PAGES* (10),VOL* (3),JOURN* (4)
*     ..
*     .. External Functions ..
      INTEGER GOTHFD
      EXTERNAL GOTHFD
*     ..
*     .. Executable Statements ..
*
* eat the first field (either C or R)
      JPOS = GOTHFD(IPOS,JUNK,JLEN,OTHERS)
*
* get the second -- this is the journal title
      JPOS = GOTHFD(JPOS+1,JOURN,JLEN,OTHERS)
*
* work out the first letter for the citation
      IF (JOURN.EQ.'cacm') THEN
          FIELD(1) = 'C'
      ELSE IF (JOURN.EQ.'toms') THEN
          FIELD(1) = 'T'
      ELSE
          FIELD(1) = 'X'
      END IF
*
* third field is pages; fourth is volume field
      JPOS = GOTHFD(JPOS+1,PAGES,PLEN,OTHERS)
      NXTOTH = GOTHFD(JPOS+1,VOL,VLEN,OTHERS)
*
* build the reference
      DO 10 I = 1,VLEN
          FIELD(I+1) = VOL(I:I)
   10 CONTINUE
      FIELD(VLEN+2) = SEPCHR
      DO 20 I = 1,PLEN
          IF (PAGES(I:I).NE.PAGSEP) THEN
              FIELD(VLEN+2+I) = PAGES(I:I)
          ELSE
              FLEN = VLEN + I + 1
              RETURN
          END IF
   20 CONTINUE
      FLEN = VLEN + PLEN + 2
      END
      SUBROUTINE OUTALG(NALG)
*
* output the complete entry for an algorithm.
* nalg gives the position within the arrays.
*
* The following latex template is used
*
* \noindent
* {\small \begin{tabular}{p{0.5in}p{2.75in}p{1.1in}}
* <algorithm label>
* &{\raggedright
* <title>
* \\}&{\raggedright {\bf
* <main reference> }
* <other references>
* \\} \end{tabular}}
*
* set up the strings and their respective lengths.
* Note: the backslash character is used extensively in these
*       strings. Its definition may be implementation dependent.
*       The following is for SUN and 4.3bsd Unix f77 compilers.
*
*     .. Parameters ..
      CHARACTER*(*) NOIND
      INTEGER LNOIND
      PARAMETER (NOIND='\\noindent',LNOIND=9)
      CHARACTER*(*) TABST
      INTEGER LTABST
      PARAMETER (TABST=
     +   '{\\small \\begin{tabular}{p{0.4in}@{}p{2.75in}@{}p{1.4in}}',
     +          LTABST=56)
      CHARACTER*(*) TABEND
      INTEGER LTABND
      PARAMETER (TABEND='\\\\}\\end{tabular}}',LTABND=17)
      CHARACTER*(*) RAGGD
      INTEGER LRAGGD
      PARAMETER (RAGGD='&{\\raggedright ',LRAGGD=15)
      CHARACTER*(*) RGIT
      INTEGER LRGIT
      PARAMETER (RGIT='\\\\}&{\\raggedright {\\bf ',LRGIT=23)
      CHARACTER BLNK,CURL
      INTEGER LBLNK,LCURL
      PARAMETER (BLNK=' ',LBLNK=1,CURL='}',LCURL=1)
      INTEGER BUFMAX,MAXLEN
      PARAMETER (BUFMAX=80,MAXLEN=160)
      INTEGER MAXALG,TITLEN,OTHLEN
      PARAMETER (MAXALG=750,TITLEN=160,OTHLEN=300)
*     ..
*     .. Scalar Arguments ..
      INTEGER NALG
*     ..
*     .. Scalars in Common ..
      INTEGER NOALGS
*     ..
*     .. Arrays in Common ..
      INTEGER LENREF(MAXALG),LENTIT(MAXALG),PAGE(MAXALG),SHSORT(MAXALG),
     +        VOL(MAXALG)
      CHARACTER JOUR(MAXALG),REFS(OTHLEN,MAXALG),TITLES(TITLEN,MAXALG),
     +          SHAR(MAXALG)*3,LABEL(MAXALG)*6
*     ..
*     .. Local Scalars ..
      INTEGER BUFLEN,FLEN
      CHARACTER BUFFER* (BUFMAX)
*     ..
*     .. Local Arrays ..
      CHARACTER ARRBUF(MAXLEN)
*     ..
*     .. External Functions ..
      INTEGER LENGTH
      EXTERNAL LENGTH
*     ..
*     .. External Subroutines ..
      EXTERNAL GETOTC,MAINRF,OUTELT,OUTLIN
*     ..
*     .. Common blocks ..
      COMMON /CALLIN/VOL,PAGE,NOALGS,LENTIT,LENREF,SHSORT
      COMMON /CALLST/TITLES,REFS,JOUR,SHAR,LABEL
*     ..
*     .. Save statement ..
      SAVE /CALLST/,/CALLIN/
*     ..
*     .. Executable Statements ..
      CALL OUTLIN(NOIND,LNOIND)
      CALL OUTLIN(TABST,LTABST)
      CALL OUTLIN(LABEL(NALG),LENGTH(LABEL(NALG)))
      CALL OUTLIN(RAGGD,LRAGGD)
      CALL OUTELT(TITLES(1,NALG),LENTIT(NALG))
      CALL OUTLIN(RGIT,LRGIT)
      CALL MAINRF(BUFFER,BUFLEN,NALG)
      BUFLEN = BUFLEN + 1
      BUFFER(BUFLEN:BUFLEN) = CURL
      CALL OUTLIN(BUFFER(1:BUFLEN),BUFLEN)
      IF (LENREF(NALG).NE.0) THEN
          CALL GETOTC(ARRBUF,FLEN,REFS(1,NALG),LENREF(NALG))
          CALL OUTELT(ARRBUF,FLEN)
      END IF
      CALL OUTLIN(TABEND,LTABND)
*
* output a blank line -- well almost
      CALL OUTLIN(BLNK,LBLNK)
      END
      SUBROUTINE OUTCUM
*
* open file for output of bibtex references
* unit number nout must not clash with nin and nerr
* which are initialized in routine setup
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
      OPEN (UNIT=NOUT,FILE='cum.tex',STATUS='unknown',FORM='formatted',
     +     ACCESS='sequential',ERR=10)
      RETURN
* deal with possible error opening output file
   10 CALL ERROR('Failed to open file for bibtex database')
      END
      SUBROUTINE OUTELT(ARRAY,LENGTH)
*
* Outputs a character array array(1..length) using a buffer
* of size at most maxbuf characters, if the length is greater
* than maxbuf splits the string on the nearest space.
*
*     .. Parameters ..
      INTEGER MAXBUF
      CHARACTER BLNK
      PARAMETER (MAXBUF=80,BLNK=' ')
*     ..
*     .. Scalar Arguments ..
      INTEGER LENGTH
*     ..
*     .. Array Arguments ..
      CHARACTER ARRAY(LENGTH)
*     ..
*     .. Local Scalars ..
      INTEGER I,J,OLDPTR,PTR
      CHARACTER BUFFER* (MAXBUF)
*     ..
*     .. External Subroutines ..
      EXTERNAL ERROR,OUTLIN
*     ..
*     .. Executable Statements ..
      OLDPTR = 1
      PTR = MAXBUF
   10 IF (PTR.GE.LENGTH) THEN
          J = 0
          DO 20 I = OLDPTR,LENGTH
              J = J + 1
              BUFFER(J:J) = ARRAY(I)
   20     CONTINUE
          CALL OUTLIN(BUFFER,LENGTH-OLDPTR+1)
          RETURN
      ELSE
          DO 30 J = PTR,OLDPTR,-1
              IF (ARRAY(J).EQ.BLNK) THEN
                  GO TO 40
              END IF
   30     CONTINUE
          CALL ERROR('Unbreakable string in outelt - increase maxbuf')
   40     PTR = J - 1
          J = 0
          DO 50 I = OLDPTR,PTR
              J = J + 1
              BUFFER(J:J) = ARRAY(I)
   50     CONTINUE
          CALL OUTLIN(BUFFER,PTR-OLDPTR+1)
          OLDPTR = PTR + 1
          PTR = PTR + MAXBUF
          GO TO 10
      END IF
      END
      SUBROUTINE OUTEND
*
* output the final lines of the LaTeX file
*     .. Parameters ..
      CHARACTER*(*) ENDLIN
      INTEGER LENDLN
      PARAMETER (ENDLIN='\\end{document}',LENDLN=14)
*     ..
*     .. External Subroutines ..
      EXTERNAL OUTLIN
*     ..
*     .. Executable Statements ..
      CALL OUTLIN(ENDLIN,LENDLN)
      END
      SUBROUTINE OUTHED
*
* output the initial LaTeX commands for the index
* document. These take the form
*
* \documentstyle{report}
* \begin{document}
* \title{Cumulative Index to the ACM Algorithms}
* \author{Tim Hopkins and David Morse\\
*         University of Kent\\
*         Canterbury\\
*         Kent, CT2 7NF, UK}
* \date{24 June 1989}
* \maketitle
*
* NOTE: This routine may contain implementation dependent
*       character definition of backslash
*
*     .. Parameters ..
      CHARACTER*(*) STYLE
      INTEGER LSTYLE
      PARAMETER (STYLE='\\documentstyle{report}',LSTYLE=22)
      CHARACTER*(*) BEGIND
      INTEGER LBEGIN
      PARAMETER (BEGIND='\\begin{document}',LBEGIN=16)
      CHARACTER*(*) TITLE
      INTEGER LTITLE
      PARAMETER (TITLE='\\title{Cumulative Index to the ACM Algorithms}'
     +          ,LTITLE=46)
      CHARACTER*(*) AUTHR1
      INTEGER LAUTH1
      PARAMETER (AUTHR1='\\author{Tim Hopkins and David Morse\\\\'//
     +          'University of Kent\\\\',LAUTH1=57)
      CHARACTER*(*) AUTHR2
      INTEGER LAUTH2
      PARAMETER (AUTHR2='Canterbury\\\\Kent, CT2 7NF, UK}',LAUTH2=30)
      CHARACTER*(*) DATE
      INTEGER LDATE
      PARAMETER (DATE='\\date{24 June 1989}',LDATE=19)
      CHARACTER*(*) MAKET
      INTEGER LMAKET
      PARAMETER (MAKET='\\maketitle',LMAKET=10)
*     ..
*     .. External Subroutines ..
      EXTERNAL OUTLIN
*     ..
*     .. Executable Statements ..
*
      CALL OUTLIN(STYLE,LSTYLE)
      CALL OUTLIN(BEGIND,LBEGIN)
      CALL OUTLIN(TITLE,LTITLE)
      CALL OUTLIN(AUTHR1,LAUTH1)
      CALL OUTLIN(AUTHR2,LAUTH2)
      CALL OUTLIN(DATE,LDATE)
      CALL OUTLIN(MAKET,LMAKET)
      END
      SUBROUTINE OUTSH(CURPT)
*
* output the classification heading -- in sharst(curpt)
* for the S section only output the heading for the main
* section. Group subsections together.
*
*
* ***** Possibly MACHINE DEPENDENT CHARACTER VARIABLE BACKSLASH
* ***** used
*
* LaTeX template used is of the form
*
* \begin{center}{\small\bf
*  <shar heading>
* }\end{center}
*
* shar heading is of the form A1:Real arithmetic
*
* For the subclassifications of the S (Special function)
* section the second line of the above template is replaced by
* just the mnenomic
*
*     .. Parameters ..
      CHARACTER*(*) CENHED
      INTEGER LHED
      PARAMETER (CENHED='\\begin{center}{\\small\\bf',LHED=24)
      CHARACTER*(*) CENEND
      INTEGER LEND
      PARAMETER (CENEND='}\\end{center}',LEND=13)
      CHARACTER*(*) SEP
      INTEGER LSEP
      PARAMETER (SEP=' : ',LSEP=3)
      INTEGER MAXBUF,MNLEN,STLEN
      PARAMETER (MNLEN=3,STLEN=45,MAXBUF=MNLEN+LSEP+STLEN)
      INTEGER NCLASS
      PARAMETER (NCLASS=57)
*     ..
*     .. Scalar Arguments ..
      INTEGER CURPT
*     ..
*     .. Scalars in Common ..
      INTEGER MAXSHR,SVAL
*     ..
*     .. Arrays in Common ..
      INTEGER SHARPT(NCLASS+1)
      CHARACTER SHARMN(NCLASS)* (MNLEN),SHARST(NCLASS)* (STLEN)
*     ..
*     .. Local Scalars ..
      INTEGER ALEN,MLEN
      CHARACTER BUFFER* (MAXBUF)
      LOGICAL NOS
*     ..
*     .. External Functions ..
      INTEGER LENGTH
      EXTERNAL LENGTH
*     ..
*     .. External Subroutines ..
      EXTERNAL OUTLIN
*     ..
*     .. Common blocks ..
      COMMON /CSHAR/SHARMN,SHARST
      COMMON /CSHARN/MAXSHR,SHARPT,SVAL
*     ..
*     .. Save statement ..
      SAVE /CSHAR/,/CSHARN/,NOS
*     ..
*     .. Data statements ..
      DATA NOS/.true./
*     ..
*     .. Executable Statements ..
      CALL OUTLIN(CENHED,LHED)
      IF (SHARMN(CURPT) (1:1).EQ.'S') THEN
*
* first time we set centred heading
          IF (NOS) THEN
              NOS = .false.
              ALEN = LENGTH(SHARST(SVAL))
              MLEN = LENGTH(SHARMN(SVAL))
              BUFFER(1:MLEN) = SHARMN(SVAL) (1:MLEN)
              BUFFER(MLEN+1:MLEN+LSEP) = SEP
              MLEN = MLEN + LSEP
              BUFFER(MLEN+1:MLEN+ALEN) = SHARST(SVAL) (1:ALEN)
              MLEN = MLEN + ALEN
              CALL OUTLIN(BUFFER(1:MLEN),MLEN)
              CALL OUTLIN(CENEND,LEND)
              CALL OUTLIN(CENHED,LHED)
          END IF
*
* otherwise set just a section number -- the point is that
* the descriptions are like books
          MLEN = LENGTH(SHARMN(CURPT))
          CALL OUTLIN(SHARMN(CURPT) (1:MLEN),MLEN)
      ELSE
*
* just set the heading
          ALEN = LENGTH(SHARST(CURPT))
          MLEN = LENGTH(SHARMN(CURPT))
          BUFFER(1:MLEN) = SHARMN(CURPT) (1:MLEN)
          BUFFER(MLEN+1:MLEN+LSEP) = SEP
          MLEN = MLEN + LSEP
          BUFFER(MLEN+1:MLEN+ALEN) = SHARST(CURPT) (1:ALEN)
          MLEN = MLEN + ALEN
          CALL OUTLIN(BUFFER(1:MLEN),MLEN)
      END IF
      CALL OUTLIN(CENEND,LEND)
      END
      SUBROUTINE SAVALG
*
* save the title and other references, shar index
* volume, journal and start page number.
* form sorted shar index - this is inefficient but
* it's not going to be done very often
*     .. Parameters ..
      INTEGER MAXALG,AUTLEN,TITLEN,KEYLEN,OTHLEN
      PARAMETER (MAXALG=750,AUTLEN=80,TITLEN=160,KEYLEN=400,OTHLEN=300)
*     ..
*     .. Scalars in Common ..
      INTEGER NOALGS,NUMBER,PAGEND,PAGEST,VOLUME,YEAR
      CHARACTER ALABEL* (6),JOURNL* (4),MONTH* (9),LANG* (3),SHARE* (3)
*     ..
*     .. Arrays in Common ..
      INTEGER LENREF(MAXALG),LENTIT(MAXALG),PAGE(MAXALG),SHSORT(MAXALG),
     +        VOL(MAXALG)
      CHARACTER AUTHOR(AUTLEN),JOUR(MAXALG),KEYWDS(KEYLEN),
     +          OTHERS(OTHLEN),REFS(OTHLEN,MAXALG),TITLE(TITLEN),
     +          TITLES(TITLEN,MAXALG),SHAR(MAXALG)*3,LABEL(MAXALG)*6
*     ..
*     .. Local Scalars ..
      INTEGER I
*     ..
*     .. External Functions ..
      INTEGER LENARR
      EXTERNAL LENARR
*     ..
*     .. Common blocks ..
      COMMON /CALLIN/VOL,PAGE,NOALGS,LENTIT,LENREF,SHSORT
      COMMON /CALLST/TITLES,REFS,JOUR,SHAR,LABEL
      COMMON /CREFNO/VOLUME,NUMBER,YEAR,PAGEST,PAGEND
      COMMON /CREFST/ALABEL,JOURNL,MONTH,LANG,SHARE,AUTHOR,TITLE,KEYWDS,
     +       OTHERS
*     ..
*     .. Save statement ..
      SAVE /CREFNO/,/CREFST/,/CALLST/,/CALLIN/
*     ..
*     .. Executable Statements ..
      NOALGS = NOALGS + 1
*
* save title and reference
      LENTIT(NOALGS) = LENARR(TITLE,TITLEN)
      LENREF(NOALGS) = LENARR(OTHERS,OTHLEN)
      DO 10 I = 1,LENTIT(NOALGS)
          TITLES(I,NOALGS) = TITLE(I)
   10 CONTINUE
      DO 20 I = 1,LENREF(NOALGS)
          REFS(I,NOALGS) = OTHERS(I)
   20 CONTINUE
*
* save volume, page and shar index
      VOL(NOALGS) = VOLUME
      PAGE(NOALGS) = PAGEST
      SHAR(NOALGS) = SHARE
      LABEL(NOALGS) = ALABEL
*
* code for journal
      IF (JOURNL.EQ.'toms') THEN
          JOUR(NOALGS) = 'T'
      ELSE IF (JOURNL.EQ.'cacm') THEN
          JOUR(NOALGS) = 'C'
      ELSE
          JOUR(NOALGS) = 'X'
      END IF
      END
      INTEGER FUNCTION SEMPOS(IPOS,OTHERS)
*
* Routine to find next semicolon in character array others
* on or after position ipos. The semicolon is used as a
* separator between references to remarks and certificates.
*
*     .. Parameters ..
      CHARACTER SEMIC
      PARAMETER (SEMIC=';')
*     ..
*     .. Scalar Arguments ..
      INTEGER IPOS
*     ..
*     .. Array Arguments ..
      CHARACTER OTHERS(*)
*     ..
*     .. Local Scalars ..
      INTEGER ITEMP
*     ..
*     .. Executable Statements ..
      ITEMP = IPOS
   10 IF (OTHERS(ITEMP).NE.SEMIC) THEN
          ITEMP = ITEMP + 1
          GO TO 10
      END IF
      SEMPOS = ITEMP
      END
      SUBROUTINE SHARFX(SHTEMP)
*
* Remove any leading or included blanks and translate all
* alphabetic character to upper case
*
*     .. Parameters ..
      CHARACTER BLNK
      PARAMETER (BLNK=' ')
*     ..
*     .. Scalar Arguments ..
      CHARACTER SHTEMP* (*)
*     ..
*     .. Local Scalars ..
      INTEGER I,IPOS,JPTR
      CHARACTER CHAR,LCSTR* (26),UCSTR* (26)
*     ..
*     .. External Functions ..
      INTEGER LENGTH
      EXTERNAL LENGTH
*     ..
*     .. External Subroutines ..
      EXTERNAL ERROR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC INDEX,LEN
*     ..
*     .. Executable Statements ..
      UCSTR = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      LCSTR = 'abcdefghijklmnopqrstuvwxyz'
*
* get rid leading or included blanks
      JPTR = 0
      DO 10 I = 1,LENGTH(SHTEMP)
          CHAR = SHTEMP(I:I)
          IF (CHAR.NE.BLNK) THEN
              JPTR = JPTR + 1
              SHTEMP(JPTR:JPTR) = CHAR
          END IF
   10 CONTINUE
*
* flag null string
      IF (JPTR.EQ.0) THEN
          CALL ERROR('Blank shar field')
      END IF
*
* convert case
      DO 20 I = 1,JPTR
          IPOS = INDEX(LCSTR,SHTEMP(I:I))
          IF (IPOS.NE.0) THEN
              SHTEMP(I:I) = UCSTR(IPOS:IPOS)
          END IF
   20 CONTINUE
*
* pad with trailing blanks
      DO 30 I = JPTR + 1,LEN(SHTEMP)
          SHTEMP(I:I) = BLNK
   30 CONTINUE
      END


