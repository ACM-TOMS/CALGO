! Time-stamp:  <2015-12-10 19:46:34 m>
!>> 2011-09-29 messy Krogh, Math a la Carte, Inc. (Copyrighted), Started code.
!>> 2014-01-16 messy Krogh Final version for publication?
!>> 2014-03-13 messy Krogh Fixed problem with $FIxx at end of text.
!>> 2014-06-26 messy Krogh Changed to print scalar Infinities with "Inf".
!>> 2014-09-14 messy Krogh Fixed off by 1 bug in sparse vector output with a 0.
!>> 2014-12-23 messy Krogh Fixed sparse vectors output when all are 0.
!>> 2015-02-07 messy Krogh Fixed small errors revealed by static analysis.
!>> 2015-06-16 messy Krogh Fixed minor bug in sparse output & $O indexes < 0
!>> 2015-12-06 messy Krogh Set up for use as include at referee's suggestion

! This code is free for noncommercial use as long as the following source files
! are included:  messy_m.f90, precision_m.f90, tmessy.f90, sample_m.f90,
! Thrdtmessy.F90, Makefile, messy_doc.tex, and messy_doc.bib.

! Used for printing error messages and for pretty output of other types.
! For documentation of use see the files mentioned above.

! Basic use: (Use _d for double, _s for single, or _q for quad precision.)

! use messy_gen, only: messy_d=>messy, rk_d=>rk ! Change _d for diff. precisions
! type(messy_ty) :: e
! Set data
! call messy(e, "Text, see first big block of comments below", ...)
! Just below "subroutine messy(...) is the list of arguments and their types.
  use iso_fortran_env, only: OUTPUT_UNIT, ERROR_UNIT
  implicit none
  integer, parameter :: numdig = ceiling(-log10(epsilon(1.0_rk)))
  type, public :: messy_ty ! Values to be set or looked at by the user.
    character (len=32) :: ename="Undefined" ! Name printed in error messages.
    integer :: fpprec = numdig ! Default for floating point precision
    integer :: kdf = numdig ! Current default real precision.
    integer :: line_len = 128 ! Default for line length
    integer :: munit = OUTPUT_UNIT ! Message unit number
    integer :: eunit = OUTPUT_UNIT ! ERROR_UNIT mixes up output with piping
    integer :: maxerr = 0 ! Max value of 1000 * (10*stop + print) + |index|
    integer :: lstop = 3 ! Stop indexes <= this don't stop
    integer :: lprint = 3 ! Print indexes <= this don't print
    integer :: errcnt = 0 ! Count of the number of error messages, incremented
!                       by 1000000 for internal errors inside messy.
    integer :: dblev = 3 ! If 0, an immediate return is made (unless text
    !                   starts with "$E"), else a $K<integer> will behave
    !                   as if reaching the end of text if <integer> is > dblev.
    integer, private :: lentry = 0 ! Tells state on first entry.
    integer, private :: ounit = OUTPUT_UNIT ! Current output unit.
    integer, private :: lenlin = 128 ! Current line length
    integer, private :: sunit = 0 ! Unit number for scratch file.
    integer, private :: tabspa = 10 ! Space fill up col.  mult. of tabspa.
    integer, private :: tbrk(7) = 0 ! At most 6 groups for table output.
    integer, private :: brk_lines = 0 ! Number of lines going to scratch file.
  end type messy_ty
contains
  recursive subroutine messy(e,text,idat,rdat,imat,rmat,zdat,zmat,ix,ptext)
    type(messy_ty), intent(inout) :: e
    character (len=*), intent(in) :: text
    character (len=*), optional,  intent(in) :: ptext
    integer, optional, intent(in) :: idat(:), imat(:,:), ix(:)
    real(rk), optional, intent(in) :: rdat(:), rmat(:,:)
    complex(rk), optional, intent(in) :: zdat(:), zmat(:,:)
! Processes Messages -- Actions are controlled by text.

! This routine is intended for use primarily by other library routines.
! Users of library routines may want to set or look at values in a variable
! declared with type messy_ty.

! The character after a '$' (or something else if parameter SC is changed)
! defines actions as follows.

!   A  Print a matrix from rmat.  If column and rows headings have not been
!      changed using $O, column headings have the form Col nnn, and Row headings
!      have the form Row nnn, where nnn is the index of the column or row.  Note
!      that the format specified by $F is never used for matrix output as it is
!      awkward to line up columns properly.
!   B  Break, ends all actions and restores all user changeable actions back
!      to their default values.
!   C  Continue, can be used to end an integer when followed by a digit.
!      When the this is the last thing in text, and we are processing an error
!      message, the error message is not ended at this point, but is continued
!      on the next call.
!   D  Used to specify a temporary number of significant digits for floating
!      point output.  Followed by characters for an integer.  Followed by a
!      non-digit restores the default and treats this character as if no $D
!      appeared.  A 0 restores the usual default value.  A negative number
!      indicates you want that number of digits after the decimal point.
!      The usual default is always restored after a $B.
!   E Start an error message.  The next two character are digits, the first
!      gives the stop index, and the next the print index.  If the stop index is
!      0, is is not treated as an error, but can be used to limit printing of
!      other messages.  If this is an error message the index of the error,
!      defined by the package generating the error message, is in {\tt ix(1)} if
!      it is present (useful if some integer vector is part of the error
!      message), and else is in {\tt idat(1)}.  (And if {\tt idat} is not
!      present, then a 0 is printed for the error index.)  Note that after any
!      print from {\tt idat} the next integer printed from {\tt idat} will come
!      from the next location.  If the stop index is 0, print will still come
!      out on the error unit.
!   F Define the alternative format for integer or floating point output.  The
!     $F is followed by a ``D'', or one of the letters: ``IFE'' followed by
!     ``digits.digits'', where the ``. digits'' is optional for ``I''.  The
!     ``D'' case is primarily intended for use with tables, see item ``H''
!     below, and gives an ``ES'' format to print the number of significant
!     digits that would ordinarily be given by e%kdf (or was last set by $D) and
!     with a width sufficient to give one space at the start of the number.  The
!     ``E'' cases are converted to Fortran's scientific edit descriptor ``ES''.
!     This also does the type of output specified by the format.  (Either the
!     next integer from idat in the ``I'' case, the next real from rdat for real
!     formats, or in the $ZF case, the next complex number from zdat.  In the
!     complex case the format provided is used for both the real and imaginary
!     part unless the format is terminated with a ``,'', in which case the
!     format for the imaginary part follows.)  These formats are not saved from
!     one call to the next.
!   G As for R below, except use the last floating format defined by a "$F",
!     above.
!   H In the middle of text, $H specifies a preferred place for breaking a line.
!     If this is before the start of text, it specifies the start of a table.
!     By a table we mean text that is arranged in columns where the caller
!     indicates what is to go in each column.  If there are more columns than
!     will fit on a line, text that does not fit will be saved in a scratch
!     file.  In most cases the end of the table will be indicated with text
!     containing nothing but the $B.  The user is responsible for setting up
!     headings, and formatting the following lines so things line up as desired.
!     The following lines will make use of formatting available through the $F
!     command above.  The heading line has the following form:
!     ``$H|nctext1|nctext2| ...  |nctextk|'', where ``|'' is any character the
!     user wants to use, different than anything else in the current text
!     string; n is a number giving the number of spaces needed for the column
!     (if 0 or missing the column is later formatted with ``$FD", see ``F''
!     above); and c is a C, R, or L indicating how the following text is to be
!     justified in the column.  The test program tmessy.F90 gives an example for
!     setting up tables.  One can have at most lenbuf-50 characters in a heading
!     where lenbuf is a parameter in messy_m.F90}, currently = 256.  Note that
!     one can make up headings, and output table entries with a size that
!     depends on the precision, by using ``$FD", and setting e%kdf to a value
!     that depends on the precision (perhaps using ``$D'').
!   I Set nidat = nidat + 1 and print idat(nidat), and continue.
!   J As for I above, except use the last integer format defined by a "$F",
!     see above.
!   K This is followed by an integer (assumed 0 with no digits given).  If that
!     integer is > e%dblev, then actions are as if text has ended at this point.
!     If e\%dblev is 0, then an immediate return is made unless text starts with
!     ``$E'', or we are already processing an error message.  The $K<integer> is
!     ignored when processing an error message.  If this feature is used, then
!     the smaller the integer following the $K, the more likely the following
!     text is to print.  The larger the number in e\%dblev, the more likely text
!     after a ``$K'' is to print.
!   L Followed by an integer gives a new line length.  Internally the number
!     specified is replaced as necessary to get it in the interval [40,
!     lenbuf-50] (parameter lenbuf is now 256).
!   M Print a matrix from imat.  Headings are as for $A.  $M$T results in the
!     transposed matrix being printed.
!   N Start a New line, and continue.
!   O Define starting indexes for vector or matrix output and for the matrix
!     case, alternative formats for output of column and row labels. The $O may
!     only appear immediately after a $A, $M, $V, or $W or after a $T
!     following one of these.  Following the $O is an optional integer giving
!     the starting index value (which may be negative).  No integer gives the
!     default, which is 1.  In the matrix case one can indicate the text to
!     output for the column and row headings, and whether and where the index
!     should be printed.  A $O following a $A or $M, is followed by text
!     which specifies the desired result for columns, a $O, then specifications
!     for rows, and then a terminating $O.  Following the optional integer
!     giving the first index, one may have nothing, or a specification of the
!     number of characters in the text for the headings in each column, and
!     where indexes are to be printed followed by the actual text used for
!     headings.  ``<#'' means print the index first with # characters in the
!     text for each column heading (all must be the same).  ``>#'' is the
!     same but the heading text follows the index.  `` #'' gives no index, but
!     just # characters of heading.  ``|##'' has the first # characters of
!     heading text then the index, and then the following # characters of
!     heading text.  The sum of the digits in the ## must not exceed 9.

!     If the text is the same for all columns (or rows) then text contains L
!     characters, otherwise it contains some multiple of L characters, and
!     different text is used for each heading until running out of text in which
!     case the last given is repeated.  One can get the default actions by not
!     using the $O, with $O$O$O, with $O>4Col $O>4Row $O, $O4Col $O$O, or
!     $O$O4Row $O.  If column headings were to be Earth, Air, Fire, and Water
!     and row labels were to have the form "Case <index>:", this would be $O
!     5Earth Air FireWater$O|51Case :$O.  Column headings are centered over the
!     data for that column, unless the column headings are wider than needed by
!     the data in which case the data is right justified with the heading.
!   P Print text in ptext.  Useful for long messages with small changes.
!   Q This is used to print integers as bit strings.  If ix is not defined it is
!     assumed that the bit string is in a single integer, and that leading 0's
!     are not to be printed.  If ix is present then $|$ix(1)$|$ defines the
!     number of bits, and if ix(1) $> 0$, then this gives the number of bits and
!     all are printed, otherwise it is assumed that leading 0's should not be
!     printed, and the bit string consists of at least one integer.  Idat
!     contains the bit strings much like when idat is used for output of
!     integers or integer vectors, with the difference that a single bit string
!     may require more than 1 word.  Following the $Q must be a B, O or a Z,
!     indicating that one wants the bits output in binary, octal, or
!     hexadecimal.  Following this must be an I or a V indicating that one wants
!     just a single bit string output, or would like to output the contents of
!     idat as a vector of bit strings.  Since Fortran does not have a bit data
!     type, the high order bit of the integers is ignored, thus if one has 32
!     bit integers, each integer holds only 31 bits.
!   R Set nrdat=nrdat+1 and print rdat(nrdat), and continue.
!   S Print the real sparse vector, with row/column indexes in idat and values
!     in rdat.  Text at the start of the line containing the $S is printed at
!     the start of each line.
!   T Tab to column that is multiple of tabspa set with $<integer>T, see 0--9
!     below.  If used immediately after a $A or $M, $T causes the transpose
!     of a matrix to be printed.
!   U Set the output unit to the following integer.  This output unit will be
!     used until the next $B.
!   V Print a floating point vector from rdat.  If other numbers have been
!     printed from rdat, then the vector starts with the first unprinted
!     number.
!   W As for $V, but for integers from idat.
!   X Print an integer from ix. $X prints ix(1), $X<integer> prints
!     ix(<integer>).
! Y,y If you would like to repeat what is in part of text ix(1) times, you can
!     do this with $Y, and then indicate the end of that block with $y. After
!     processing the text in $Y...$y ix(1) times, text is then processed in the
!     usual way.  This was introduced primarily for tables ($H) so that tables
!     could have a variable number of columns.  When used for column headings,
!     the $Y and the $y should both be followed by the special character used to
!     delimit column headings.
!   Z  Used for complex output.  Should be followed immediately by R,
!      F, G, V, or A, where these letters cause actions just like these
!      ordinarily do for the real case.
! 0-9  Starts <integer> and then either a F, G, J, T, or a ' '. will repeat the
!      $F, $G or $J action that number of times, or set the spacing for tabs or
!       output that many blanks.
! <sc> a single <sc> (currently a '$') is printed, continue till the next <sc>.
!      else given an error.

!
! ************************** Internal Variables ************************
!
! atext  Formal argument to arrchk giving text for array being checked.  First
!   letter gives the type, and rest gives the array name.
! bigfrac Fractional part (base 10) of rmax.
! bnd    Formal argument to get_row_label giving bounds for current array.
! bndi   Upper bound for idat.
! bndm1  Upper bound on first dimension for a matrix.
! bndm2  Upper bound on second dimension for a matrix.
! bndoth Used in process_matrix for bounds on rows/columns.
! bndr   As for bndi, but for rdat.
! bndx   As for bndi, but for ix.
! bndz   As for bndi, but for zdat.
! e%brk_lines Number of lines going to scratch file.
! buf    Character string holding characters to be output.
! bz     Number of bits in a single word of idat.
! c      Used for temp. storage of a character
! c1     Temporary character when processing complex scalar
! cs     Character selected for main action.
! csh    Character used to delimit headings
! dlen   Number of bits used to construct a single output character.
! e      Argument giving the messy_ty structure in the call.
! errchar Parameter giving the error message. separator.  Default is a '$'.
!   '0' give no separator.
! errdat Used to hold array dimensions for an error message.
! e%ename Package name to print with error messages.

! e%errcnt Count of error messages.
! errtxt Text passed as part of internal error in this code.
! etext  Text used as template for internal error messages to messy.
! e%eunit  Unit number for output of error messages.
! fmtt   Temporary for format for both integer and floating point output.
! fmtg   Format set by user for floating point output.
! fmtz   Format set by user for complex output.
! fmtj   Format set by user for integer output.
! e%fpprec Default value for kdf.
! got0   Set true if we have seen a zero value and no negatives.
! got0i  As for got0 for the imaginary parts of complex numbers.
! hj     Character used to indicate justification of headings.
! hlen   Length of the current heading
! hmap   Array used in getting character to print for hex digits.
! i      Temporary index
! iachar0 Parameter giving the value of iachar('0')
! iadjust Holds index adjustment for columns and rows.
! ib     Used to hold bits needed for first word of bit strings.
! ib1    Used for bits before the output of the first group.
! icol   Index tracking columns in a matrix.
! id     Number of decimal digits to print after the decimal point.
! idat   Formal array, containing integer data to output.
! ierr   Formal argument to mess_error giving the index associated with an
!   internal error.
! ifmt   Index tracking where things are stored in fmtt.
! igroup Number of output characters between blanks in binary/hex output.
! imag   Magnitude of integer to output, with negative sign if integer
!   is < 0.
! imat   Formal array, containing integer 2D array to output.
! e%inerr  0 if not processing an error message, 1 if printing an error
!   message, -1 if in an error message that is not being printed, and >1
!   if printing an error message that stops.  Set to -2 when the error
!   message is supposed to stop.
! iout   Integer to be output.
! irc    = 1 for rows, = 2 for columns when determining labels for
!   matrix output.
! itb    Index into tbrk to track where (of if) things are on the scratch file.
! itemp  Array used to store integers about to be output.  (Used if got0=.true.)
! iw     Width to be used in a floating pt. or integer format.
! iw1    Saves value of iw for real part of complex number.
! iwg    Value like iw for user set real format.
! iwj    Value like iw for user set integer format.
! iwmin  Minimal width for an output number due to column headings.
! iwz    Value like iw for user set complex format.
! ix     Formal array, containing alternative integer data to output.
! j      Used as a temporary index.
! k      Used as a temporary index.
! k1     Used as a temporary index.
! kbexp  Base 10 exponent of rmax.
! e%kdf  Current number of digits to print for floating point.  If <= 0, then
!   -e%kdf gives the number of digits to print after the decimal point.
! keb    Used in tracking number of characters needed to represent the base
!   10 exponent for rmax.
! kf1    Last character index used in krfmt(1)
! klen   Temp. value for dlen on output of first digit.
! krxxx   Names starting with kr are all arrays with two elements.  The first
!         is for columns, the second is for rows.
! kr1    Characters before the numeric label.  < 0 if no label.
! kr2    Digits after the numeric label.
! krfmt  Formats for column and row headings.
! krm    Start of next text for column and row labels
! krt1   Where first text for column/row labels is stored in krfmt.    
! krt2   Where last text for column/row labels is stored in krfmt.    
! krw    Number of character used for column and row labels.
! ks     Used +/- 1 depending on sign in getint
! e%kscrn  Number of lines to "print" before pausing.
! ksexp  Base 10 exponent of rmin.
! kskip  Number of leading words for bit output that are all 0 bits.
! l      Used as a temporary index.
! lbuf   Position of characters in BUF, usually the last to print.
! lbufh  Position in buffer just prior to heading text.
! lenbuf Parameter giving the number of character in BUF.
! e%lenlin Gives number of character in output lines.
! lentxt Length of character input array text.
! lentxth Value of lentxt when lentxt is set to location of next special
!   character use in headings.
! e%lentry Tells what to do on entry (and sometimes other places.)
!   = 0  Has not been called yet.
!   = 1  Processing usual messages.
!   = 2  Processing an error message that does not stop.
!   = 3  Processing an error message that stops
!   = 4  Set just prior to stopping on an error message.
!   = 5  Set for an error message that does not print.
!   =-1  Working on later lines for a table.
!   =-2  Processing a table heading, saw first $H.
!   =-3  Saw the second $H for a table heading
! linebrk A desired value of lbuf for breaking a line.
! e%line_len Gives the default line length.
! e%lprint For error messages with a print level <= e%lprint nothing is
!   printed (unless the message would result in a stop).
! e%lstop  As for lprint, except with the stop level, and stopping.
! e%maxerr Maximum value seen for 1000 * (10*stop + print) + index.
! maxv   Maximum value in an integer vector.
! e%munit  Output unit used for messages that aren't in an error message.
! n      Used for number of numbers to print on line.
! nbits  Number of bits needed to represent bit/hex data.
! neg    Set to 1 if a there is a negative number
! negi   Same as neg, but for the imaginary part of complex numbers.
! nidat  Index of last item in idat printed.
! nrdat  Index of last item in rdat printed.
! ntext  text(ntext+1:?) contains the next text data.
! ntextr text(ntextr+1:?) contains place to start over on a repeat.
! num    Number of times some type of user formatted output is to be repeated.
! e%numdig Global default real precision
! numr   Number of times actions from text are to be repeated.
! nwords Number of integers (words) used to hold bit/hex data.
! nzdat  Index of last item in zdat printed.
! e%ounit Index of the current output unit.
! p      Used for stop index in err_head.
! ptext   Formal argument for extra text to output.
! rdat   Formal array, containing real data to output.
! rmat   Optional formal array containing a matrix.
! rmax   Largest | value | among real numbers to be printed.
! rmaxi  Same as rmax, but for the imaginary part of complex numbers.
! rmin   Smallest nonzero ! value | for among real numbers to be printed.
! rmini  Same as rmin, but for the imaginary part of complex numbers.
! rndup  Constant array used in testing if numbers will require more digits
!   if rounded up for the current number of digits output..
! rout   Value of a real number to print.
! rtemp  Array used to store real numbers about to be output, as for itemp.
! s      Used for stop index in err_head.
! sc     Parameter for special character used to introduce actions.
!   Default value is '$'.  Changing sc does not change the $'s printed as part
!   of error messages.  These can be changed by changing errchar.
!   (Note that using "\" for sc is not portable.)
! smafrac Fractional part (base 10) of rmin.
! e%sunit  Index for the scratch unit, 0 if not yet assigned.
! e%tabspa $T tabs to next col. that is multiple of tabspa, always gives at
!   least one blank.
! e%tbrk Stores column indexes where line breaks occur when processing tables.
!   tbrk(1:2) are also used for the stop and print indexes for error messages.
! text   Formal argument giving the character string from which all text
!   is taken.
! trans  Set to .true. if matrix is to be transposed, else .false.
! vlabsz Size of space needed for a row label.
! zdat   Formal array, containing complex data to output.
! zmat   Formal array, containing complex 2D array to output.
!
! ************************** Variable Declarations *********************
!
    integer, parameter :: lenbuf=250
! You can change the $ in the line below, but if you do you should change $
! to the same character wherever a $ appears here except when the $ is 
! contained in the string "!$OMP" (there are very few of these).
    character(len=1), parameter :: sc = '$'
    character(len=1), parameter :: errchar = '$'
    integer, parameter :: iachar0=iachar('0'), iachar9=iachar('9')
    character (len=77) :: etext="$E48Error in inputs to messy.  The offending&
     & input text is near the end of:$C"
    integer :: errdat, lentxt, lentxth, nidat, nrdat, ntextr, nzdat
    integer, volatile :: ntext ! Volatile to get around bug in unnamed compiler
    character (len=13) :: fmtg ! User set formats
    character (len=40) :: fmtz
    character (len=8) :: fmtj
    integer :: iwg, iwj, iwz ! iw for user set formats
    character (len=40) :: fmtt
    character (len=64) :: krfmt(2)
    character :: c, c1, cs, csh, hj

    character (len=lenbuf) :: buf ! Buffer used in constructing lines.
    integer :: bndi, bndm1, bndm2, bndoth, bndr, bndx, bndz, hlen, i, id, ifmt,&
      & iout, itb, iw, iw1, iwmin, j, k, k1, kbexp, keb, ks, ksexp, l, lbuf,&
      & lbufh, linebrk, maxv, n, iadjust(2), neg, negi, num, numr, p, s,&
      & vlabsz, itemp(80)

    integer :: icol, kf1, kr1(2), kr2(2), krm(2), krt1(2), krt2(2), krw(2)

    real(rk) :: bigfrac, rmax, rmaxi,  rmin, rmini, rout, rtemp(50), smafrac
    logical :: got0, got0i, trans
    real(rk), parameter :: rndup(numdig) =&
      & [ (1.0_rk - .5_rk/10._rk**k, k=1, numdig) ]
    integer :: bz, dlen, ib, ib1, igroup, klen, kskip, nbits, nwords

#if numt_ > 1
! A little space is saved by leaving numt_ undefined if you never use openmp.
    save bigfrac,bndi,bndoth,bndm1, bndm2,bndr,bndx,bndz,buf,bz,&
     & c,c1,cs,csh,dlen,errdat,fmtt,fmtg,fmtj,fmtz,got0,got0i,trans,&
     & hj,hlen,i,iadjust,ib,ib1,icol,id,ifmt,igroup,iout, itb,&
     & itemp,iw,iw1,iwg,iwj,iwmin,iwz,j,k,k1,kbexp,keb,kf1,klen,kr1,kr2,&
     & krfmt,krm,krt1,krt2,krw,ks,ksexp,kskip,l,lbuf,lbufh,lentxt,&
     & lentxth,linebrk,maxv,n,nbits,neg,negi,nidat,nrdat,ntext,ntextr,num,&
     & numr,nwords,nzdat,p,rmax,rmaxi,rmin,rmini,rout,rtemp,s,&
     & smafrac,vlabsz
!$OMP THREADPRIVATE(bigfrac,bndi,bndoth,bndm1,bndm2,bndr,bndx,bndz,buf,bz,&
!$OMP c,c1,cs,csh,dlen,errdat,fmtt,fmtg,fmtj,fmtz,got0,got0i,trans,&
!$OMP hj,hlen,i,iadjust,ib,ib1,icol,id,ifmt,igroup,iout, itb,&
!$OMP itemp,iw,iw1,iwg,iwj,iwmin,iwz,j,k,k1,kbexp,keb,kf1,klen,kr1,kr2,&
!$OMP krfmt,krm,krt1,krt2,krw,ks,ksexp,kskip,l,lbuf,lbufh,lentxt,&
!$OMP lentxth,linebrk,maxv,n,nbits,neg,negi,nidat,nrdat,ntext,ntextr,num,&
!$OMP numr,nwords,nzdat,p,rmax,rmaxi,rmin,rmini,rout,rtemp,s,&
!$OMP smafrac,vlabsz)
#endif

!
! ************************* Start of Executable Code *******************
!
!
    if (e%dblev == 0) then
      if (e%lentry <= 1) then
        if (index(text, "$E") == 0) return
      end if
    end if
    if (e%lentry <= 0) then
      if (e%lentry < 0) then ! Processing a table.
        itb = 2 ! Index to tbrk for finding line breaks.
        e%lentry = -1
        lbufh = 0
        e%lenlin = e%tbrk(2)
        e%ounit = e%munit
      else ! Very first call, set some defaults.
        e%ounit = e%munit
        e%lenlin = e%line_len
        e%lentry = 1
        e%kdf = e%fpprec
      end if
    else if (e%lentry == 5) then ! Non-printing error message.
      if (index(text, "$C") == 0) then ! No more continuations
        e%ounit = e%munit
        e%kdf = e%fpprec
      end if
      e%lentry = 1
      return
    end if
    lentxt = len(text) ! Initializations on entry
    fmtt(:) = ' ' ! Initialization of the entire arrays or the NAG
    fmtj(:) = ' ' ! compiler nags about it.
    fmtg(:) = ' '
    fmtz(:) = ' '
    linebrk = 100000
    ntext = 0
    nidat = 0
    nrdat = 0
    nzdat = 0
    bndi = -1000001
    bndr = -1000001
    bndx = -1000001
    bndz = -1000001
    iwmin = 0
    lbuf = 0
    numr = 0
    outer: do ! For going from one '$' to the next
      do while (lbuf > e%lenlin) ! Flush what we can
        call bufout
      end do ! End of flushing all that could output a line.
      if (ntext >= lentxt) then ! End of input
        if (e%lentry < -1) then
          if (buf(lbuf:lbuf) == csh) then
            lbuf = lbuf - 1
          end if
          call col_head
        end if
      end if
      if (ntext >= lentxt) then
        call done
        return
      end if
      k = index(text(ntext+1:lentxt), sc)
      if (k /= 1) then
        if (k == 0) then
          k = min(lenbuf-lbuf, lentxt - ntext)
          buf(lbuf+1:lbuf+k) = text(ntext+1:ntext+k)
          lbuf = lbuf + k
          ntext = ntext + k
          cycle
        end if
        k = min (k-1, lentxt - ntext, lenbuf - lbuf)
        buf(lbuf+1:lbuf+k) = text(ntext+1:ntext+k)
        ntext = ntext + k
        lbuf = lbuf + k
        if (lbuf > e%lenlin) cycle
      end if
      ntext = ntext + 2
      if (ntext > lentxt) then
        call mess_error(4, " ")
        cycle
      end if
      cs = text(ntext:ntext)
      num = 1
      select case (cs) ! Start of all the different cases
      case ('A') ! Print a real matrix
        if (.not. present(rmat)) then
          call mess_error(2, "Armat")
          return
        end if
        bndm1 = ubound(rmat,1)
        bndm2 = ubound(rmat,2)
        rmax = 0.0_rk
        rmin = huge(1.0_rk)
        neg = 0
        got0 = .false.
        do j = 1, bndm2
          do i = 1, bndm1
            if (rmat(i,j) < 0.0_rk) then
              rmax = max(rmax, -rmat(i,j))
              rmin = min(rmin, -rmat(i,j))
              neg = 1
            else if (rmat(i,j) > 0.0_rk) then
              rmax = max(rmax, rmat(i,j))
              rmin = min(rmin, rmat(i,j))
            else if (rmat(i,j) == 0.0_rk) then
              got0 = .true.
            end if ! NaN'a will take care or themselves
          end do
        end do
        if (neg > 0) got0 = .false.
        call process_matrix
        lbuf = 0

      case ('B') ! Break, end of text line and action
        if (e%lentry < 0) call finish_table
        call done
        e%lentry = 1
        e%ounit = e%munit
        e%kdf = e%fpprec
        e%lenlin = e%line_len
        return

      case ('C') ! Just continue (special if an error message)
        if (ntext == lentxt) then
          if (e%lentry > 1) then
            if (lbuf > 0) then
              write (e%eunit, '(A)') buf(1:lbuf)
            end if
            return
          end if
        end if

      case ('D') ! Specify temporary digits for floating point output.
        e%kdf = min(getint(), numdig)
        if ((e%kdf == 0) .and. (text(ntext:ntext) /= '0')) e%kdf = e%fpprec

      case ('E') ! Start an error message
        if (e%lentry < 0) call finish_table
        if (lbuf > 0) then
          write(e%ounit, '(A)') buf(1:lbuf)
          lbuf = 0
          call mess_error(1, "$E must be the first thing in text.")
          return
        end if
        ntext = ntext + 2
        if (ntext > lentxt) then
          call mess_error(4, " ")
          return
        end if
        e%tbrk(1) = iachar(text(ntext-1:ntext-1)) - iachar0
        e%tbrk(2) = max(e%tbrk(1), iachar(text(ntext:ntext)) - iachar0)

        if ((e%tbrk(1) < 0) .or. (e%tbrk(1) > 9) .or.&
          & (e%tbrk(2) > 9)) then
          call mess_error(1, "Two digits must follow a $E")
          return
        end if
        if (present(ix)) then
          e%tbrk(3) = ix(1)
        else
          e%tbrk(3) = 0
          if (e%tbrk(1) > 0) then
            if (present(idat)) then
              if (arrchk("Iidat")) then
                nidat = nidat + 1
                if (nidat <= bndi) e%tbrk(3) = idat(nidat)
              end if
            end if
          end if
        end if
        if (e%tbrk(1) <= e%lstop) e%maxerr = max(e%maxerr,&
          & 1000* (10 * e%tbrk(1) + e%tbrk(2)) + abs(e%tbrk(3)))
        if (e%tbrk(2) <= e%lprint) then
          if (index(text, "$C") /= 0) e%lentry = 5
          return
        end if
        e%lentry = 2
        if (e%tbrk(1) > e%lstop) e%lentry = 3
        if (e%tbrk(1) /= 0) then
          write(e%eunit, '()') ! Output blank line.
          e%ounit = e%eunit
          call err_head
          i=1
          do while (iachar(e%ename(i:i)) > iachar(' '))
            buf(i:i) = e%ename(i:i)
            i = i + 1
          end do
          write(e%eunit,'(A," reports error: Stop level=",i0,&
            &" Print level=",i0," Error index=", i0)')&
            & buf(1:i-1), e%tbrk(1), e%tbrk(2), e%tbrk(3)
        end if
      case ('F') ! Alternative format for integer or floating point output.
        call f_output
        if (ntext > lentxt) return

      case ('G') ! Print rdat(nrdat+1) using format specified with 'F' above.
        call g_output
        if (ntext > lentxt) return

      case ('H') ! Process heading for table line breaks.
        if (e%lentry > 0) then ! Process a table
          if (lbuf > 0) then ! Setting line break location
            linebrk = lbuf
          else
            if (ntext >= lentxt-5) return
            e%lentry = -2 ! Process text to get heading and line breaks.
            itb = 1
            ntext = ntext + 1
            csh = text(ntext:ntext)
            lentxth = lentxt
            lbufh = lbuf
            call col_head
          end if

        else if (e%lentry == -1) then
          call mess_error(1, &
            "Are you starting a table without ending a previous one with $B?")
          return
        end if
      case ('I') ! Print idat(nidat)
        if (bndi == -1000001) then
          if (.not. arrchk("Iidat")) return
        end if
        iw = 1
        neg = 0
        nidat = nidat + 1
        if (nidat > bndi) then
          errdat = bndi
          call mess_error(5, "idat")
          return
        end if
        iout = idat(nidat)
        if (iout < 0) then
          neg = 1
        else if (iout == 0) then
          lbuf = lbuf + 1
          buf(lbuf:lbuf) = '0'
          cycle outer
        end if
        maxv = abs(iout)
        call get_formati
        write (buf(lbuf+1:lbuf+iw), fmtt) iout
        lbuf = lbuf + iw

      case ('J') ! Print idat(nidat using format defined by 'F' above
        call j_output
        if (ntext > lentxt) return

      case ('K') ! Stop if following integer is < e%dblev
        if (getint() > e%dblev) then
          if (e%lentry <= 1) lentxt = ntext
        end if

      case ('L') ! Set the current line length
        k1 = getint()
        e%lenlin = min(lenbuf-40, max(40, k1))

      case ('M') ! Print an integer matrix
        if (.not. present(imat)) then
          call mess_error(2, "Mimat")
          return
        end if
        bndm1 = ubound(imat,1)
        bndm2 = ubound(imat,2)
        iw = 2
        neg = 0
        maxv = 0
        got0 = .false.
        do j = 1, bndm2
          do i = 1, bndm1
            if (imat(i,j) < 0) then
              neg = 1
              maxv = max(maxv, -i,imat(i,j))
            else if (imat(i,j) > 0) then
              maxv = max(maxv, imat(i,j))
            else ! Must be 0, but don't want a -0.
              got0 = .true.
            end if
          end do
        end do
        if (neg > 0) got0 = .false.
        call process_matrix
        lbuf = 0

      case ('N') ! Start a new line
        write(e%ounit, '(A)') buf(1:lbuf)
        linebrk = 100000
        lbuf = 0

      case ('O') ! Alternative labels for matrix output.
        call mess_error(1,&
          & '$O is only allowed as part of vector or matrix output.')
        return

      case ('P') ! Output text in ptext
        if (.not. present(ptext)) then
          call mess_error(2, "Pptext")
          return
        end if
        k = min (len(ptext), lenbuf-lbuf)
        buf(lbuf+1:lbuf+k) = ptext
        lbuf = lbuf + k

      case ('Q') ! Binary or Hex output ($Q{B,H,O}{I,V}
        ntext = ntext + 2
        if (ntext > lentxt) then
          call mess_error(4, " ")
          return
        end if
        if (bndi == -1000001) then
          if (.not. arrchk("Iidat")) return
        end if
        bz = bit_size(idat(1)) - 1
        if (present(ix)) then
          nbits = ix(1)
          if (nbits == 0) nbits = -bz
          nwords = 1 + (abs(nbits) -1) / bz
        else
          nwords = 1
          nbits = -bit_size(idat(nidat+1)) + 1
        end if
        c = text(ntext-1:ntext-1)
        select case(c) ! Test for valid input
        case ('B')
          dlen = 1
        case ('O')
          dlen = 3
        case ('Z')
          dlen = 4
        case default
          call mess_error(11,"Letter following $Q must a B, O or Z.&
            &  You have a ")
          return
        end select
        c1 = text(ntext:ntext)
        kskip = 0
        if (nbits > 0) then
          ib = nbits - bz * (nwords-1)
        else
          if (c1 /= 'V') then ! Output of a single bit map
            do
              ib = idat(nidat+kskip+1)
              call bit_count(ib)
              kskip = kskip + 1
              if ((ib>0) .or. (nwords<=kskip)) exit
            end do
            ib = min(ib, -nbits - bz * (nwords-kskip))
          else ! A vector of bit maps
            do
              ib = 0 ! Use this and do loop below if no Fortran 2008
              do i = nidat+kskip+1, bndi, nwords
                ib = ior(ib, idat(i))
              end do ! Instead of this do, use line below if Fortran 2008.
!                ib = iany(idat(nidat+kskip+1:bndi:nwords)) ! (Fortran 2008)
              call bit_count(ib)
              kskip = kskip + 1
              if ((ib>0) .or. (nwords<=kskip)) exit
            end do
          end if
          kskip = kskip - 1
          ib = max(ib, 1)
        end if
        if (nidat+nwords > bndi) then
          errdat = bndi
          call mess_error(5, "idat")
          return
        end if
        nwords = nwords - kskip
        nbits = ib + bz * (nwords-1)
        iw = 1 + (ib + bz*(nwords-1)-1)/dlen ! Space needed for the bit string
        igroup = iw
        ib1 = ib
        if (iw >= 10) then ! Extra space for blank separators
          iw = iw + (iw-1)/8
          igroup = 8
          ib1 = 1 + mod(nbits-1, 8*dlen)
        end if
        iw = iw + 3 ! For the {BZ}" at start and " at the end
        select case (c1)
        case ('I')
          if (c1 == 'I') then ! A single bit map
            k = linebrk
            if (k > lbuf) k = -2
            vlabsz = mod(lbuf-k+2+ib1, igroup+1) - 1
            if (lbuf > e%lenlin) call bufout
          end if
          call bits_out ! Output of one bit/hex map.
        case ('V')! A vector of bit maps
          iw = iw + 1 ! For the extra blank between vector entries.
          krfmt(2)(1:12) = '(SS,i1,": ")' ! Get format for row labels
          vlabsz = 3
          l = nidat + 1
          if (l < 0) then
            krfmt(2)(3:3) = 'P'
            vlabsz = 4
          end if
          n = (bndi - nidat) / nwords
          k = 10
          do while (k < n)
            k = k * 10
            vlabsz = vlabsz + 1
          end do
          write (krfmt(2)(6:6), '(ss,i1)') vlabsz - 2
          if (n * nwords /= (bndi - nidat)) then
            call mess_error(9, "")
          end if ! Keep going and print what we can.
          if (n <= 0) then
            buf(lbuf+1:lbuf+16) = "Vector is empty."
            n = 0
          else if (lbuf == 0) then
            n = 0
          else ! n is number of bit strings, iw space for a bit string.
            if (n * iw > (e%lenlin - lbuf)) then ! Need > 1 line
              if (lbuf + 2*iw <= e%lenlin) then ! >0 bit maps will fit.
                k = vlabsz
                do while (k < lbuf) ! Find place to start print
                  k = k + iw !       So that all lines up
                end do
                if (lbuf < k) then
                  if (k < e%lenlin) then ! Extra spaces so things line up
                    buf(lbuf+1:k) = ' '
                    lbuf = k
                  end if
                end if
              else if (lbuf /= 0) then
                write (e%ounit, '(A)') buf(1:lbuf)
                lbuf = 0
              end if
            end if
          end if
          do
            if (lbuf == 0) then
              write (buf(1:vlabsz), krfmt(2)(1:12)) 1+(nidat+1-l)/nwords
              lbuf = vlabsz
            end if
            call bits_out
            if (nidat > bndi - nwords) cycle outer
            ib = nbits - bz * (nwords-1)
            if (lbuf + iw < e%lenlin) then
              lbuf = lbuf + 1
              buf(lbuf:lbuf) = ' '
              cycle
            end if
            write (e%ounit, '(A)') buf(1:lbuf)
            lbuf = 0
          end do
        case default
          c = c1
          call mess_error(11, "Second letter following $Q must be an I or a W.&
            &  You have a ")
        end select

      case ('R') ! Print rdat(nrdat+1) with default format
        if (bndr == -1000001) then
          if (.not. arrchk("Rrdat")) return
        end if
        nrdat = nrdat + 1
        if (nrdat > bndr) then
          errdat = bndr
          call mess_error(5, "rdat")
          return
        end if
        rmax = rdat(nrdat)
        if (rmax < 0.0_rk) then
          rmax = -rmax
          neg = 1
        else if (rmax > 0.0_rk) then
          neg = 0
        else if (rmax == 0.0_rk) then ! Just print a 0.0 for a 0.
          lbuf = lbuf + 2
          buf(lbuf-1:lbuf) = '0.'
          cycle outer
        else ! Assume it is a Nan
          buf(lbuf+1:lbuf+3) = "NaN"
          lbuf = lbuf + 3
          cycle outer
        end if
        if (rmax > huge(1.0_rk)) then
          if (neg /= 0) then
            lbuf = lbuf + 1
            buf(lbuf:lbuf) = '-'
          end if
          buf(lbuf+1:lbuf+3) = "Inf"
          lbuf = lbuf + 3
          cycle outer
        end if
        rmin = rmax
        ifmt = 1
        fmtt(1:1) = '('
        call get_format
        if (lbuf+iw >= lenbuf) call bufout
        write (buf(lbuf+1:lbuf+iw), fmtt) rdat(nrdat)
        lbuf = lbuf + iw
      case ('S') ! Output real sparse vector
        if (bndi == -1000001) then
          if (.not. arrchk("Iidat")) return
        end if
        if (bndr == -1000001) then
          if (.not. arrchk("Rrdat")) return
        end if
        if (bndi /= bndr) then
          call mess_error(8, " ")
          return
        end if
        maxv = maxval(idat(1:bndi))
        iw = 1
        call get_formati
        fmtt(1:13) = "(A,99(:' (',i"
        iw1 = iw
        write(fmtt(14:15), "(ss,i2)") iw
        fmtt(16:20) = ",',',"
        ifmt = 20
        call vec_setup
        fmtt(ifmt:ifmt+5) = "')'))"
        iw = iw + iw1 + 3
        n = (e%lenlin - lbuf)/ iw
        if (n == 0) then ! To avoid looping
          call bufout
          n = e%lenlin/ iw
        end if
        do
          k = min(n, bndi-nidat, bndr-nrdat)
          if (k == 0) exit
          if (got0) then
            rtemp(1:k) = abs(rdat(nrdat+1:nrdat+k))
            write (e%ounit, fmtt) buf(1:lbuf), (idat(i+nidat), rtemp(i), i=1, k)
          else
            write (e%ounit, fmtt) buf(1:lbuf), (idat(i+nidat), rdat(i+nrdat),&
              & i=1, k)
          end if
          nrdat = nrdat + k
          nidat = nidat + k
        end do
        lbuf = 0
        linebrk = 100000
      case ('T') ! Tab
        k = e%tabspa - mod(lbuf, e%tabspa)
        if (k == 0) k = e%tabspa
        buf(lbuf+1:lbuf+k) = ' '
        lbuf = lbuf + k

      case ('U') ! Set new ounit
        e%ounit = getint()
        if (e%ounit == 0) e%ounit = e%munit

      case ('V') ! Print a floating point vector from rdat
        fmtt(1:8) = "(99(:1x,"
        ifmt = 8
        call vec_setup
        call get_row_label(bndr)
        if (lbuf > e%lenlin+16) then
          call bufout
        end if
        n = bndr - nrdat
        if (n <= 0) then
          buf(lbuf+1:lbuf+16) = "Vector is empty."
          n = 0
        else if (lbuf == 0) then
          n = 0
        end if
        if (n * iw > (e%lenlin - lbuf)) then ! Shift so things will line up
          k = vlabsz ! The vector is going to require more than one line.
          do while (k < lbuf)
            k = k + iw
          end do
          if (lbuf < k) then ! Extra spaces so things line up.
            buf(lbuf+1:k) = ' '
            lbuf = k
          end if
          n = (e%lenlin - lbuf) / iw
        end if
        if (n > 0) then
          if (got0) then
            rtemp(1:n) = abs(rdat(nrdat+1:nrdat+n))
            write (buf(lbuf+1:lbuf+n*iw), fmtt) rtemp(1:n)
          else
            write (buf(lbuf+1:lbuf+n*iw), fmtt) rdat(nrdat+1:nrdat+n)
          end if
        end if
        if (lbuf+n > 0) write (e%ounit, '(A)') buf(1:lbuf+n*iw)
        lbuf = 0
        nrdat = nrdat + n
        if (nrdat >= bndr) cycle outer
        krfmt(2)(15:) = fmtt(2:)
        n = (e%lenlin - vlabsz) / iw
        do
          k = min(nrdat+n, bndr)
          if (got0) then
            rtemp(1:k-nrdat) = abs(rdat(nrdat+1:k))
            write (e%ounit, krfmt(2)) nrdat+iadjust(2), rtemp(1:k-nrdat)
          else
            write (e%ounit, krfmt(2)) nrdat+iadjust(2), rdat(nrdat+1:k)
          end if
          nrdat = k
          if (nrdat == bndr) cycle outer
        end do
        linebrk = 100000

      case ('W') ! Print an integer vector from idat
        if (bndi == -1000001) then
          if (.not. arrchk("Iidat")) return
        end if
        iw = 2
        neg = 0
        maxv = 0
        got0 = .false.
        do i = nidat+1, bndi
          if (idat(i) < 0) then
            neg = 1
            maxv = max(maxv, -idat(i))
          else if (idat(i) > 0) then
            maxv = max(maxv, idat(i))
          else ! Must be 0, but don't want a -0.
            got0 = .true.
          end if
        end do
        if (neg > 0) got0 = .false.
        call get_formati
        call get_row_label(bndi)
        if (lbuf > e%lenlin+16) then
          call bufout
        end if
        n = bndi - nidat
        if (n <= 0) then
          buf(lbuf+1:lbuf+16) = "Vector is empty."
          n = 0
        else if (lbuf == 0) then
          n = 0
        end if
        if (n * iw > (e%lenlin - lbuf)) then ! Shift so things will line up
          k = vlabsz ! The vector is going to require more than one line.
          do while (k < lbuf) ! Get start at least to lbuf.
            k = k + iw
          end do
          if (lbuf < k) then ! Extra spaces so things line up
            buf(lbuf+1:k) = ' '
            i = abs(idat(1)) ! Check if first number would fit in the space
            j = k-1
            if (idat(1) < 0) j = j - 1
            do while (i >= 10) 
              i = i / 10
              j = j - 1
            end do
            if (j > lbuf) then ! We have space to fit this number in.
              nidat = nidat + 1
              write(buf(j+1:k), '(i0)') idat(nidat)
            end if
          end if
          n = (e%lenlin - k) / iw
          lbuf = k
        end if
        if (n > 0) then
          if (got0) then
            itemp(1:n) = abs(idat(nidat+1:nidat+n))
            write (buf(lbuf+1:lbuf+n*iw), fmtt) itemp(1:n)
          else
            write (buf(lbuf+1:lbuf+n*iw), fmtt) idat(nidat+1:nidat+n)
          end if
        end if
        if (lbuf+n > 0) write (e%ounit, '(A)') buf(1:lbuf+n*iw)
        lbuf = 0
        nidat = nidat + n
        if (nidat >= bndi) cycle outer
        krfmt(2)(15:) = fmtt(2:)
        n = (e%lenlin - vlabsz) / iw
        do
          k = min(nidat+n, bndi)
          if (got0) then
            itemp(1:k-nidat) = abs(idat(nidat+1:k))
            write (e%ounit, krfmt(2)) nidat+iadjust(2), itemp(1:k-nidat)
          else
            write (e%ounit, krfmt(2)) nidat+iadjust(2), idat(nidat+1:k)
          end if
          nidat = k
          if (nidat == bndi) cycle outer
        end do
        linebrk = 100000

      case ('X') ! Output integer from ix().
        if (bndx == -1000001) then
          if (.not. arrchk("Xix")) return
        end if
        k = 1
        if (ntext <= lentxt) then
          k = max(getint(), 1)
        end if
        if (k > bndx) then
          errdat = bndx
          call mess_error(5, "ix")
          return
        end if
        maxv = abs(ix(k))
        neg = 0
        if (ix(k) < 0) neg = 1
        iw = 1
        call get_formati
        write (buf(lbuf+1:lbuf+iw), fmtt) ix(k)
        lbuf = lbuf + iw
      case ('Y') ! Repeat actions from text ix(1) times.
        if (bndx == -1000001) then
          if (.not. arrchk("Xix")) return
        end if
        numr = ix(1)
        ntextr = ntext
        if (e%lentry < -1) then
          if (text(ntextr+1:ntextr+1) == csh) then
            ntext = ntext + 1
            ntextr = ntext
            call col_head
          end if
        end if
      case ('y') ! End of repeat actions
        numr = numr - 1
        if (numr > 0) then
          ntext = ntextr
        else
          if (e%lentry < -1) ntext = ntext + 1
        end if
        if (e%lentry < -1) call col_head

      case ('Z') ! For output of complex data
        ntext = ntext + 1
        c = text(ntext:ntext)
        do
          select case (c)
          case('R') ! Single complex number
            if (bndz == -1000001) then
              if (.not. arrchk("Zzdat")) return
            end if
            lbuf = lbuf + 1
            buf(lbuf:lbuf) = '('
            nzdat = nzdat + 1
            if (nzdat > bndz) then
              errdat = bndz
              call mess_error(5, "zdat")
              return
            end if
            rout = real(zdat(nzdat)) ! First we deal with the real part
            c1 = ','
            do
              rmax = abs(rout)
              rmin = rmax
              neg = 0
              if (rout == 0.0_rk) then
                buf(lbuf+1:lbuf+2) = "0."
                lbuf = lbuf + 2
              else if (rmax >= rmin) then ! Not a NaN
                if (rout < 0) neg = 1
                fmtt(1:1) = '('
                ifmt = 1
                call get_format
                if (lbuf+iw >= lenbuf) then
                  call bufout
                end if
                write (buf(lbuf+1:lbuf+iw), fmtt) rout
                lbuf = lbuf + iw
              else
                buf(lbuf+1:lbuf+3) = "NaN"
                lbuf = lbuf + 3
              end if
              lbuf = lbuf + 1
              buf(lbuf:lbuf) = c1
              if (c1 .eq. ')') cycle outer
              c1 = ')'
              rout = aimag(zdat(nzdat)) ! Next the imaginary part
            end do
            cycle outer
          case ('F') ! Specify user format for complex output.
            fmtz(1:8) = "(99('(',"
            j = 9
            kf1 = ntext
            iwz = 3
            do l = 1, 2
              if (c /= ',') ntext = kf1
              ntext = ntext + 1
              c = text(ntext:ntext)
              if ((c == 'F') .or. (c == 'E') .or. (c == 'D')) then
                call get_user_format(fmtz)
                iwz = iwz + k1
              else
                call mess_error(11, "Format letter for $ZF must be 'F', 'E',&
                  & or 'D'.  You have a ")
                return
              end if
              if (l == 2) exit
              fmtz(j:j+4) = ",',',"
              j = j + 5
              ntext = ntext + 1
              c = text(ntext:ntext)
            end do
            fmtz(j:j+6) = ",')'))"
            c = 'G'
            cycle
          case ('G') ! Single complex number, user format
            if (bndz == -1000001) then
              if (.not. arrchk("Zzdat")) return
            end if
            if (fmtz(1:1) /= '(') then
              call mess_error(7, "$ZG (or F)")
              return
            end if
            nzdat = nzdat + 1
            if (nzdat > bndz) then
              call mess_error(5, "zdat")
              return
            end if
            if (lbuf + iwz > e%lenlin) then
              call bufout
            end if

            write(buf(lbuf+1:lbuf+iwz+1), fmtz) zdat(nzdat:nzdat)
            lbuf = lbuf + iwz
            cycle outer
          case ('V') ! Output a complex vector from zdat
            if (bndz == -1000001) then
              if (.not. arrchk("Zzdat")) return
            end if
            iw = -1
            rmax = 0.0_rk
            rmin = huge(1.0_rk)
            rmaxi = 0.0_rk
            rmini = huge(1.0_rk)
            neg = 0
            negi = 0
            got0 = .false.
            got0i = .false.
!!              fmtt
            do i = nzdat+1, bndz
              if (real(zdat(i)) < 0.0_rk) then
                rmax = max(rmax, -real(zdat(i)))
                rmin = min(rmin, -real(zdat(i)))
                neg = 1
              else if (real(zdat(i)) > 0.0_rk) then
                rmax = max(rmax, real(zdat(i)))
                rmin = min(rmin, real(zdat(i)))
              else
                got0 = .true.
              end if
              if (aimag(zdat(i)) < 0.0_rk) then
                rmaxi = max(rmaxi, -aimag(zdat(i)))
                rmini = min(rmini, -aimag(zdat(i)))
                negi = 1
              else if (aimag(zdat(i)) > 0.0_rk) then
                rmaxi = max(rmaxi, aimag(zdat(i)))
                rmini = min(rmini, aimag(zdat(i)))
              else
                got0i = .true.
              end if
            end do
            if (neg > 0) got0 = .false.
            if (negi > 0) got0i = .false.
            fmtt(1:12) = '(99(:1x,"(",'
            ifmt = 12
            call get_format
            fmtt(ifmt:ifmt+4) = '")")))'
            iw = iw + iw1 + 4
            call get_row_label(bndz)
            if (lbuf > e%lenlin+16) then
              call bufout
            end if
            n = bndz - nzdat
            if (n <= 0) then
              buf(lbuf+1:lbuf+16) = "Vector is empty."
              n = 0
            else if (lbuf == 0) then
              n = 0
            end if
            if (n * iw > (e%lenlin - lbuf)) then ! Shift so things will line up
              k = vlabsz ! The vector is going to require more than one line.
              do while (k < lbuf)
                k = k + iw
              end do
              if (lbuf < k) then
                buf(lbuf+1:k) = ' '
                lbuf = k
                n = (e%lenlin - lbuf) / iw
              end if
              if (n > 0) then
                if (got0) then
                  rtemp(1:2*n:2) = abs(real(zdat(1:n)))
                else
                  rtemp(1:2*n:2) = real(zdat(1:n))
                end if
                if (got0i) then
                  rtemp(2:2*n:2) = abs(aimag(zdat(1:n)))
                else
                  rtemp(2:2*n:2) = aimag(zdat(1:n))
                end if
                write (buf(lbuf+1:lbuf+n*iw), fmtt) rtemp(1:2*n)
              end if
            end if
            if (lbuf+n > 0) write (e%ounit, '(A)') buf(1:lbuf+n*iw)
            lbuf = 0
            nzdat = nzdat + n
            if (nzdat >= bndz) cycle outer
            krfmt(2)(15:) = fmtt(2:)
            n = (e%lenlin - vlabsz) / iw
            if (n == 0) then
              neg = iw + vlabsz ! lenlin must be this big
              call mess_error(10, "")
              return
            end if
            do
              k = min(nzdat+n, bndz)
              k1 = 2 * (k - nzdat)
              if (got0) then
                rtemp(1:k1) = abs(real(zdat(nzdat+1:k)))
              else
                rtemp(1:k1:2) = real(zdat(nzdat+1:k))
              end if
              if (got0i) then
                rtemp(2:k1:2) = abs(aimag(zdat(nzdat+1:k)))
              else
                rtemp(2:k1:2) = aimag(zdat(nzdat+1:k))
              end if
              write (e%ounit, krfmt(2)) nzdat+1, rtemp(1:k1)
              nzdat = k
              if (nzdat == bndz) cycle outer
            end do
            linebrk = 100000
            cycle outer
          case ('A') ! Output a complex matrix
            if (.not. present(zmat)) then
              call mess_error(2, "Zzmat")
              return
            end if
            bndm1 = ubound(zmat,1)
            bndm2 = ubound(zmat,2)
            rmax = 0.0_rk
            rmin = huge(1.0_rk)
            rmaxi = 0.0_rk
            rmini = huge(1.0_rk)
            neg = 0
            negi = 0
            got0 = .false.
            do j = 1, bndm2
              do i = 1, bndm1
                if (real(zmat(i,j)) < 0.0_rk) then
                  rmax = max(rmax, -real(zmat(i,j)))
                  rmin = min(rmin, -real(zmat(i,j)))
                  neg = 1
                else if (real(zmat(i,j)) > 0.0_rk) then
                  rmax = max(rmax, real(zmat(i,j)))
                  rmin = min(rmin, real(zmat(i,j)))
                else if (real(zmat(i,j)) == 0.0_rk) then
                  got0 = .true.
                end if ! NaN'a will take care or themselves
                if (aimag(zmat(i,j)) < 0.0_rk) then
                  rmaxi = max(rmaxi, -aimag(zmat(i,j)))
                  rmini = min(rmini, -aimag(zmat(i,j)))
                  negi = 1
                else if (aimag(zmat(i,j)) > 0.0_rk) then
                  rmaxi = max(rmaxi, aimag(zmat(i,j)))
                  rmini = min(rmini, aimag(zmat(i,j)))
                else if (aimag(zmat(i,j)) == 0.0_rk) then
                  got0i = .true.
                end if ! NaN'a will take care or themselves
              end do
            end do
            if (neg > 0) got0 = .false.
            if (negi > 0) got0i = .false.
            call process_matrix
            lbuf = 0
            cycle outer
          case default ! Error
            call mess_error(11, "An R, F, G, V, or A, must follow $Z.&
              &  You have a ")
          end select
        end do
      case ('0':'9') ! Multiple outputs or tab settings
        ntext = ntext - 1
        num = getint()
        ntext = ntext + 1
        cs = text(ntext:ntext)
        select case (cs)
        case ('J')
          call j_output
        case ('G')
          call g_output
        case ('F')
          call f_output
        case ('T')
          e%tabspa = num
        case (' ')
          buf(lbuf+1:min(e%lenlin,lbuf+num)) = ' '
          lbuf = lbuf + num
        case default
          c = cs
          call mess_error(11, "J, G, F, T or ' ' must follow a $X.&
            &  You have a ")
          return
        end select

      case (sc) ! Just output the sc (special character.
        lbuf = lbuf + 1
        buf(lbuf:lbuf) = sc

      case default ! Output the sc followed by the character.
        c = cs
        call mess_error(11, "Following a $, must be A-Z, y, 0-9 or a $.&
          &  You have a ")
        return

      end select
    end do outer ! For going from one '$' to the next
    return
  contains

    integer function getint()
      ks = 1
      getint = 0
      if (ntext >= lentxt) return
      if (text(ntext+1:ntext+1) == '-') then
        ks = -1
        ntext = ntext + 1
      end if
      do ntext = ntext + 1, lentxt
        id = iachar(text(ntext:ntext)) - iachar0
        if ((id < 0) .or. (id > 9)) exit
        getint = 10*getint + id
      end do
      ntext = ntext - 1
      getint = getint * ks
      return
    end function getint


    subroutine bufout
      if (e%lentry >= 0) then
        l = lbuf
        if (l > linebrk) l = linebrk
        if (l > e%lenlin) then
          l = e%lenlin
          do j = e%lenlin, e%lenlin/2, -1 ! Find a place to break the line.
            if (buf(j+1:j+1) == ' ') then
              l = j
              exit
            end if
          end do
        end if
        write (e%ounit, '(A)') buf(1:l)
        do l = l+1, lbuf
          if (buf(l:l) /= ' ') exit
        end do
        buf(1:lbuf-l+1) = buf(l:lbuf) ! Shift data in buf left
        lbuf = lbuf - l + 1
        linebrk = 10000
      else
        if (itb == 1) then
          call mess_error(1, "Two columns must appear before the end of line")
          return
        else
          if (e%tbrk(itb) == 0) then
            call mess_error(1, "$H column is too wide for the line_length")
          end if
          write (e%ounit, '(A)') buf(1:e%tbrk(itb))
          if (e%tbrk(itb) < lbuf) then ! If test needed due to bug in
            buf(e%tbrk(1)+1:e%tbrk(1)+lbuf-e%tbrk(itb))=buf(e%tbrk(itb)+1:lbuf)
          end if ! NAG compiler  (Code inside the if block always needed)
          lbufh = e%tbrk(1) + (lbufh - e%tbrk(itb))
          lbuf = e%tbrk(1) + (lbuf - e%tbrk(itb))
          
          itb = itb + 1
          if (itb > 7) then
            call mess_error(1, "$H table is way too long")
            return
          end if
          if (e%tbrk(itb) /= 0) then
            if (e%lentry == -1) e%lenlin = e%tbrk(itb)
          else
            if (e%lentry == -1) then ! Scratch file has the start of line
              lbuf = 0
            end if
          end if
          if (e%ounit /= e%sunit) then
            if (e%sunit == 0) then
! If you get compilation errors on the next statement replace using s%unit=99
              open(status="SCRATCH", NEWUNIT=e%sunit, iostat=k)
!              s%sunit = 99 ! and uncomment these two lines
!              open(status="SCRATCH", unit=s%sunit, iostat=k)
! If you are using threads, you need s%sunit to depend on omp_get_thread_num().
              if (k /= 0) then
                call mess_error(1, "Error on open of needed scratch file.")
                return
              end if
            end if
            e%ounit = e%sunit
          end if
        end if
      end if
      return
    end subroutine bufout

    subroutine err_head
      k = 0
      if ((errchar /= '0') .and. (e%tbrk(1) > 0)) then
        buf(1:e%lenlin) = repeat(errchar,e%lenlin) ! Error mess. separator
        k = e%lenlin
      end if
      if (e%lentry == 4) then
        if (e%maxerr == 0) then
          write(e%eunit, '("There have been no errors prior to this one.")')
        else
          s = e%maxerr / 10000
          e%maxerr = e%maxerr - 10000 * s
          p = e%maxerr / 1000
          i = e%maxerr - 1000 * p
          write(e%eunit, '("Previously there have been ", i0, " error messages.&
            &  The most serious had a"/ "stop index of ", i0,",&
            & a print index of ", i0, ", and an |error index| of ", i0)')&
            & e%errcnt, s, p, i
        end if


        buf(5:37)= ' Fatal error -- Program stopped. '
        if (errchar == '0') then
          buf(1:4) = ' '
          k = 37
        end if
        write(e%eunit, '(A)') buf(1:k)
        stop
      end if
      if (k /= 0) write(e%eunit, '(A)') buf(1:k)
      return
    end subroutine err_head

    subroutine get_row_label(bnd)
      integer, intent(in) :: bnd ! first index - 1, and last index
      vlabsz = 2
      krfmt(2)(1:14) = '(SS,i?,":",SS,'
      iadjust(2) = 1
      if (ntext+3 <= lentxt) then
        if (text(ntext+1:ntext+2) == "$O") then
          ntext = ntext + 2
          k = getint()
          if ((k == 0) .and. (text(ntext:ntext) /= '0')) k = 1
          iadjust(2) = k ! we add 1 less than next index
        end if
      end if
      j = max(abs(iadjust(2)), bnd+iadjust(2)-1)
      k = 10
      do while (k <= j)
        k = k * 10
        vlabsz = vlabsz + 1
      end do
      if (iadjust(2) < 0) then
        vlabsz = vlabsz + 1
        krfmt(2)(3:3) = 'P'
      end if
      write (krfmt(2)(6:6), '(ss,i1)') vlabsz - 1 ! (i?:)
      return
    end subroutine get_row_label

    subroutine get_formati
      do while (maxv > 9)
        maxv = maxv / 10 ! Divides as multiplies are tricky near huge.
        iw = iw + 1
      end do
      iw = iw + neg
      iw = max(iw, iwmin)
      write(fmtt, '("(99i", i0, ")")') iw
      return
    end subroutine get_formati        

    subroutine get_format
      do
        if (rmax < rmin) then ! Have NaN's, could be big numbers
          rmax = huge(1.0_rk)
          rmin = huge(1.0_rk)
        end if
        rmax = min(rmax, huge(1.0_rk)) ! Avoiding problems with infinity
        rmin = min(rmin, rmax)
        bigfrac = log10(rmax)
        if (bigfrac >= 0.0_rk) then
          kbexp = int(bigfrac)
        else
          kbexp = int(bigfrac) - 1
        end if
        keb = 0
        if (e%kdf <= 0) then ! Fixed number of digits after the decimal.
          id = -e%kdf
          bigfrac = bigfrac + .500000001 * 10.0_rk**e%kdf
          iw = max(iwmin, id + max(0, int(bigfrac)+1) + 1 + neg)
          write (fmtt(ifmt+1:ifmt+7), '("F", i2, ".", i2, ")")') iw, id
          ifmt = ifmt + 7
          if (ifmt /= 19) return
        else
          bigfrac = bigfrac - kbexp
          if (bigfrac > rndup(e%kdf)) keb = 1
          smafrac = log10(rmin)
          if (smafrac >= 0.0_rk) then
            ksexp = int(smafrac)
          else
            ksexp = int(smafrac) - 1
            smafrac = smafrac - ksexp
          end if
          if ((kbexp < 5) .and. (ksexp > -3)) then ! I don't like big
            id = max(0, e%kdf - ksexp - 1)
            iw = max(e%kdf, kbexp+1) + kbexp - ksexp + 1 - min(ksexp+1, 0)
            if (iw <= e%kdf + 3 + keb) then ! Use an 'F' format
              iw = max(iwmin, iw + neg + keb)
              write (fmtt(ifmt+1:ifmt+7), '("F", i2, ".", i2, ")")') iw, id
              ifmt = ifmt + 7
              if (ifmt /= 19) return ! Returns unless on first complex
              go to 100 ! Avoiding labeled blocks till Fortran catches up
            end if
          end if
          iw = max(iwmin, neg + e%kdf + 4)
          kbexp = kbexp + keb
          if (max(abs(ksexp), abs(kbexp)) < 10) then
            write (fmtt(ifmt+1:ifmt+10), '("es", i2, ".", i2, "e1)")') iw, e&
              &%kdf - 1
            ifmt = ifmt + 10
            if (ifmt /= 22) return ! Returns unless on first complex
          else
            iw = iw + 1
            write (fmtt(ifmt+1:ifmt+8), '("es", i2, ".", i2, ")")') iw, e%kdf-1
            ifmt = ifmt + 8
            if (ifmt /= 20) return ! Returns unless on first complex
          end if
        end if
100     continue
        fmtt(ifmt:ifmt+4) = ',",",'
        ifmt = ifmt + 4
        rmax = rmaxi
        rmin = rmini
        neg = negi
        iw1 = iw
      end do
    end subroutine get_format

    subroutine get_user_format(vfmt)
      character(len=*) :: vfmt

      if (c /= 'D') then
        vfmt(j:j) = c
        j = j + 1
        if (c == 'E') then
          vfmt(j:j) = 'S'
          j = j + 1
        end if
        k = ntext
        k1 = getint()
        if (k1 <= 0) then
          call mess_error(1, "Field width for a format statement must be > 0")
          k1 = 10
          return
        end if
        ntext = ntext + 1
        if (ntext > lentxt) then
          if (c == 'I') return ! $FIxx is the last thing in text
          go to 10
        end if
        if (text(ntext:ntext) == '.') then
          do ntext=ntext+1,lentxt ! Using iachar for weird char. orders.
            if ((iachar(text(ntext:ntext)) < iachar0) .or.&
              & (iachar(text(ntext:ntext)) > iachar9)) exit
          end do
        else if (c /= 'I') then
          go to 10
        end if
        ntext = ntext - 1
        vfmt(j:j+ntext-k+1) = text(k+1:ntext)
        j = j + ntext - k + 1
        if (vfmt(j-1:j-1) == '.') then
          vfmt(j:j) = '0' ! Want at least one digit after the '.'.
          j = j + 1
        end if
        vfmt(j:j) = ')'
        return
      else
        k1 = e%kdf + 8
        write (fmtg(j:), '("es",I0,".",I0,")")') k1, e%kdf
        j = len_trim(fmtg)
        return
      end if
10    continue
      call mess_error(1, "Problem with text of user format.")
      return
    end subroutine get_user_format

    subroutine f_output
      j = 4
      ntext = ntext + 1
      c = text(ntext:ntext)
      select case(c)
      case ('I')
        fmtj(1:3) = "(99"
        call get_user_format(fmtj)
        iwj = k1
        call j_output
        return
      case ('F', 'E', 'D')
        fmtg(1:3) = "(99"
        call get_user_format(fmtg)
        iwg = k1
      case default ! Flag an error!  (Also error if k1 too small or big?)
        call mess_error(11, "Letter following $F must be 'I', 'F', 'E' or 'D'.&
          &  You have a ")
        return
      end select
      call g_output
      return
    end subroutine f_output

    subroutine g_output
      if (fmtg(1:1) /= '(') then
        call mess_error(7, "G")
        return
      end if
      if (bndr == -1000001) then
        if (.not. arrchk("Rrdat")) return
      end if
      if (nrdat + num > bndr) then
        errdat = bndr
        call mess_error(5, "rdat")
        return
      end if
      do i = 1, num
        if (lbuf+iwg > e%lenlin) call bufout
        nrdat = nrdat + 1
        rmax = rdat(nrdat)
        if (rmax == 0.0_rk) rmax = 0.0_rk ! Avoid printing -0's.
        write(buf(lbuf+1:lbuf+iwg), fmtg) rmax
        lbuf = lbuf + iwg
      end do
    end subroutine g_output

    subroutine j_output
      if (fmtj(1:1) /= '(') then
        call mess_error(7, "J")
        return
      end if
      if (bndi == -1000001) then
        if (.not. arrchk("Iidat")) return
      end if
      if (nidat + num > bndi) then
        errdat = bndi
        call mess_error(5, "idat")
        return
      end if
      do i = 1, num
        if (lbuf+iwj > e%lenlin) call bufout
        nidat = nidat + 1
        write(buf(lbuf+1:lbuf+iwj), fmtj) idat(nidat:nidat)
        lbuf = lbuf + iwj
      end do
    end subroutine j_output

    subroutine vec_setup
      if (bndr == -1000001) then
        if (.not. arrchk("Rrdat")) return
      end if
      rmax = 0.0_rk
      rmin = huge(1.0_rk)
      neg = 0
      got0 = .false.
      do i = nrdat+1, bndr
        if (rdat(i) == 0.0_rk) then
          got0 = .true.
        else if (rdat(i) < 0.0_rk) then
          rmax = max(rmax, -rdat(i))
          rmin = min(rmin, -rdat(i))
          neg = 1
        else if (rdat(i) > 0.0_rk) then
          rmax = max(rmax, rdat(i))
          rmin = min(rmin, rdat(i))
        end if ! NaN'a will take care or themselves
      end do
      if (neg > 0) got0 = .false.
      if (rmax == 0.0_rk) then ! Format when all is 0
        fmtt(ifmt+1:ifmt+6) = 'f3.0))'
        ifmt = ifmt+5
        iw = 4
      else
        call get_format
        fmtt(ifmt+1:ifmt+1) = ')'
        iw = iw + 1
      end if
    end subroutine vec_setup

    subroutine process_matrix ! Get info. for possibly user defined labels.
      trans = .false.
      if (ntext+2 <= lentxt) then
        if (text(ntext+1:ntext+2) == "$T") then
          ntext = ntext + 2
          trans = .true.
          i = bndm1
          bndm1 = bndm2
          bndm2 = i
        end if
      end if
      bndoth = bndm2
      do i = 1, 2 ! 1 for columns and 2 for rows
        krm(i) = 0 ! Defaults
        kr1(i) = 4
        kr2(i) = 0
        iadjust(i) = 0
        if (ntext+4 <= lentxt) then
          if (text(ntext+1:ntext+2) == "$O") then
            ntext = ntext + 2
            select case (text(ntext+1:ntext+1))
            case ('-', '0':'9')
              iadjust(i) = getint() - 1
            end select
            k = iachar(text(ntext+2:ntext+2)) - iachar0
            krm(i) = ntext + 2
            select case (text(ntext+1:ntext+1))
            case ('<')
              kr1(i) = 0
              kr2(i) = k
            case ('>')
              kr1(i) = k
              kr2(i) = 0
            case ('|')
              ntext = ntext + 1
              krm(i) = ntext + 2
              kr1(i) = k
              kr2(i) =  iachar(text(ntext+2:ntext+2)) -  iachar0
              k = kr1(1) + kr2(1)
              if (k > 9) then
                call mess_error(1, "Headings can have at most 9 characters")
                return
              end if
            case (' ')
              kr1(i) = -1
              kr2(i) = k
            case ('$')
              krm(i) = 0
            case default
              call mess_error(1, "Confused by text following $O")
              return
            end select
            k = index(text(ntext:lentxt), "$O")
            if (k == 0) then
              call mess_error(1, "$O missing in matrix format")
              return
            end if
            ntext = ntext + k
            if (i .eq. 1) ntext = ntext - 2
          end if
        end if
        krfmt(i)(1:4) = '(SS,'
        krw(i) = 0
        if (kr1(i) >= 0) krw(i) = krw(i) + 1
        k = bndoth + iadjust(i)
        do while (k >= 10)
          k = k / 10
          krw(i) = krw(i) + 1
        end do
        if (1 + iadjust(i) < 0) then
          krw(i) = krw(i) + 1
          krfmt(i)(3:3) = 'P'
        end if
        k = 4 ! To track last character in krfmt
        if (i == 1) then! The column case needs a bit more work
          iwmin = krw(1) + max(0, kr1(1)) + kr2(1) + 1
          if (cs == 'M') then ! Integer case
            call get_formati
          else
            if (rmax < rmin) then ! All 0 or NaN's
              rmax = 5000.0_rk
              rmin = 0.0_rk
            end if
            if (cs == 'A') then ! Real case
              fmtt(1:7) = "(99(1x,"
              ifmt = 7
              call get_format
              fmtt(ifmt+1:ifmt+1) = ')'
              iw = iw + 1
            else !              complex case
              fmtt(1:12) = '(99(:1x,"(",'
              ifmt = 12
              call get_format
              fmtt(ifmt:ifmt+4) = '")")))'
              iw = iw + iw1 + 4
            end if
          end if
! Set up formats
! krfmt(1) = (S?,[??x][,i?][,"text"][,i?][,"text"][,??x])
! krfmt(2) = (S?,[i?,]["text",][i?,]["text"]),SS
          j = (iw - iwmin + 1) / 2
          if (j > 0) then
            write (krfmt(1)(5:10), '(i2,"x,")') j
            k = 10
          end if
        end if
        if (kr1(i) >= 0) then
          if (kr1(i) == 0) then ! Index comes out first
            krfmt(i)(k+1:k+3) = 'i1,' ! This 1 gets reset later
            k = k + 3
          else
            if (krm(i) == 0) then ! The default case
              if (i == 1) then
                krfmt(1)(k+1:k+10) = '"Col ",i1,' ! This 1 also reset
              else
                krfmt(2)(k+1:k+10) = '"Row ",i1,' ! This 1 also reset
              end if
              k = k + 10
            else
              krfmt(i)(k+1:k+1) = '"'
              krt1(i) = k + 1 ! Where the first text is stored
              k = k + kr1(i) + 6
              krfmt(i)(k-4:k) = '",i1,' ! And this 1 etc.
            end if
          end if
          krfmt(i)(k-1:k-1) = achar(iachar0+krw(i)) ! The resetting
        end if
        if (kr2(i) > 0) then
          krfmt(i)(k+1:k+1) = '"'
          krt2(i) = k + 1 ! Where the second text is stored
          k = k + kr2(i) + 3
          krfmt(i)(k-1:k) = '",'
        end if
        if (i .eq. 1) then ! More special for the column format
          j = (iw - iwmin) / 2
          if (j > 0) then
            write (krfmt(i)(k+1:k+4), '(i2,"x,")') j
            k = k + 4
          end if
          krfmt(1)(k:k) = ')' ! Done with setting up the column format
          kf1 = k
        else ! Special for the row format
          krfmt(2)(k+1:k+3) = "SS,"
          krfmt(2)(k+4:) = fmtt(2:) ! Done with setting up the row label format
        end if
        krw(i)=krw(i)+max(kr1(i),0)+kr2(i) ! Fix final width
        bndoth = bndm1
      end do
      if (lbuf > krw(2)) then
        call bufout
      end if
      n = (e%lenlin - krw(2)) / iw
      if (n == 0) then
        neg = iw + krw(2) ! lenlin must be this big
        call mess_error(10, "")
        return
      end if
      do j = 1, bndm2, n! Output the matrix
        buf(lbuf+1:krw(2)+1) = '' ! Set up for the column headings
        lbuf = krw(2) + 1
        l = min(j + n - 1, bndm2)
        do icol = j, l
          if (krm(1) > 0) then
            if (kr1(1) > 0) then
              if (krm(1) + kr1(1) > lentxt) go to 999 ! 999 is an error exit.
              krfmt(1)(krt1(1)+1:krt1(1)+kr1(1)) = text(krm(1)+1:krm(1)+kr1(1))
              krm(1) = krm(1) + kr1(1)
            end if
            if (kr2(1) > 0) then
              if (krm(1) + kr2(1) > lentxt) go to 999
              krfmt(1)(krt2(1)+1:krt2(1)+kr2(1)) = text(krm(1)+1:krm(1)+kr2(1))
              krm(1) = krm(1) + kr2(1)
            end if
            if (krm(1) + 2 > lentxt) go to 999
            if (text(krm(1)+1:krm(1)+2) == "$O") krm(1) = 0
          end if
          if (kr1(1) < 0) then
            write(buf(lbuf+1:lbuf+iw), krfmt(1)(1:kf1))
          else
            write(buf(lbuf+1:lbuf+iw), krfmt(1)(1:kf1)) icol + iadjust(1)
          end if
          lbuf = lbuf + iw
        end do
        write (e%ounit, '(A)') buf(1:lbuf)
        lbuf = krw(2) ! Column labels are out
        buf(1:lbuf) = ' '
        do i = 1, bndm1
          if (krm(2) > 0) then
            if (kr1(2) > 0) then
              if (krm(2) + kr1(2) > lentxt) go to 999
              krfmt(2)(krt1(2)+1:krt1(2)+kr1(2)) = text(krm(2)+1:krm(2)+kr1(2))
              krm(2) = krm(2) + kr1(2)
            end if
            if (kr2(2) > 0) then
              if (krm(2) + kr2(2) > lentxt) go to 999
              krfmt(2)(krt2(2)+1:krt2(2)+kr2(2)) = text(krm(2)+1:krm(2)+kr2(2))
              krm(2) = krm(2) + kr2(2)
            end if
            if (krm(2) + 2 > lentxt) go to 999
            if (text(krm(2)+1:krm(2)+2) == "$O") krm(2) = 0
          end if
          if (cs == 'M') then ! The integer case
            if (trans) then
              itemp(1:l-j+1) = imat(j:l,i)
            else
              itemp(1:l-j+1) = imat(i,j:l)
            end if
            if (got0) itemp(1:l-j+1) = abs(itemp(1:l-j+1))
            if (kr1(2) < 0) then
              write(e%ounit, krfmt(2)) itemp(1:l-j+1)
            else
              write(e%ounit, krfmt(2)) i + iadjust(2), itemp(1:l-j+1)
            end if
          else if (cs == 'A') then ! The real case
            if (trans) then
              rtemp(1:l-j+1) = rmat(j:l,i)
            else
              rtemp(1:l-j+1) = rmat(i,j:l)
            end if
            if (got0) rtemp(1:l-j+1) = abs(rtemp(1:l-j+1))
            if (kr1(2) < 0) then
              write(e%ounit, krfmt(2)) rtemp(1:l-j+1)
            else
              write(e%ounit, krfmt(2)) i + iadjust(2), rtemp(1:l-j+1)
            end if
          else ! The complex case
            k1 = 2 * (l - j + 1)
            if (trans) then
              rtemp(1:k1:2) = real(zmat(j:l,i))
              rtemp(2:k1:2) = aimag(zmat(j:l,i))
            else
              rtemp(1:k1:2) = real(zmat(i,j:l),rk)
              rtemp(2:k1:2) = aimag(zmat(i,j:l))
            end if
            if (got0) rtemp(1:k1:2) = abs(rtemp(1:k1:2))
            if (got0i) rtemp(2:k1:2) = abs(rtemp(2:k1:2))
            if (kr1(2) < 0) then
              write(e%ounit, krfmt(2)) rtemp(1:k1)
            else
              write(e%ounit, krfmt(2)) i + iadjust(2), rtemp(1:k1)
            end if
          end if
        end do
      end do
      linebrk = 100000
      iwmin = 0
      return
999   continue ! Get here on error on runaway ntext
      call mess_error(1,&
        & "Attempt to access text past the end.  Bad count on heading?")
      return
    end subroutine process_matrix

    logical function arrchk(atext)
      character (len=*), intent(in) :: atext
      select case (atext(1:1))
      case ('I')
        arrchk = present(idat)
        if (arrchk) then
          bndi = ubound(idat,1)
          nidat = 0
          if (nidat < bndi) return
          errdat = bndi
          return
        end if
      case ('R')
        arrchk = present(rdat)
        if (arrchk) then
          bndr = ubound(rdat,1)
          nrdat = 0
          if (nrdat < bndr) return
          errdat = bndr
          return
        end if
      case ('Z')
        arrchk = present(zdat)
        if (arrchk) then
          bndz = ubound(zdat,1)
          nzdat = 0
          if (nzdat < bndz) return
          errdat = bndz
          return
        end if
      case ('X')
        arrchk = present(ix)
        if (arrchk) then
          bndx = ubound(ix,1)
          return
        end if
      end select
      call mess_error(2,atext)
      return
    end function arrchk

    subroutine bits_out ! Output of bit strings
      character(len=16), parameter :: hmap = "0123456789ABCDEF"
      k1 = 1 + (ib1 - 1) / dlen ! Number of digits to print.
      k = nidat + kskip
      nidat = k + 1
      lbuf = lbuf + 2
      buf(lbuf-1:lbuf-1) = c
      buf(lbuf:lbuf) = '"'
      if (iw+vlabsz > e%lenlin) then ! Long strings need more attention
        j = mod(k1,8)
        if (j /= 0) then ! Extra space so things line up.
          buf(lbuf+1:lbuf+8-j) = ' '
          lbuf = lbuf + 8 - j
        end if
      end if
      klen = dlen
      if (mod(nbits,klen) /= 0) klen = mod(nbits,klen)
      do
        if (klen - ib > 0) then ! Awkward case (bits in different integers)
          j = klen - ib
          i = ishft(ibits(idat(nidat),0,klen-j), j)
          nidat = nidat + 1
          ib = bz - j
          i = i + ibits(idat(nidat), ib, j)
        else
          ib = ib - klen
          i = ibits(idat(nidat),ib, klen)
        end if
        lbuf = lbuf + 1
        buf(lbuf:lbuf) = hmap(i+1:i+1) ! Output next digit
        klen = dlen ! After first use the standard length
        if (ib <= 0) then
          if (nidat -k >= nwords) exit
          nidat = nidat + 1
          ib = bz + klen - 1
        end if
        k1 = k1 - 1
        if (k1 == 0) then
          if (lbuf + 8 >= e%lenlin) then
            if (linebrk < lbuf) then
              call bufout
            else
              buf(lbuf+1:lbuf+2) = ' &'
              write(e%ounit, "(A)") buf(1:lbuf+2)
              linebrk = 100000
              lbuf = vlabsz+1
              buf(1:1) = '&'
              buf(2:lbuf) = " "
            end if
          end if
          lbuf = lbuf + 1
          buf(lbuf:lbuf) = ' '
          k1 = igroup
        end if
      end do
      lbuf = lbuf + 1
      buf(lbuf:lbuf) = '"'
      return
    end subroutine bits_out

    subroutine bit_count(ib) ! Counts number of bits to print.
      integer, intent(inout) :: ib
      k = ib
      ib = 0
      do
        if (k == 0) exit
        k = ishft(k,-1)
        ib = ib + 1
      end do
      ib = min(ib, bz)
      return
    end subroutine bit_count


    subroutine done
      if (lbuf /= 0) then
        write(e%ounit, '(A)') buf(1:lbuf)
      end if
      if (e%lentry /= 1) then
        if (e%lentry > 1) then ! Finished an error message
          e%lentry = e%lentry + 1
          call err_head
          e%lentry = 1
          write(e%eunit, '()')
          e%errcnt = e%errcnt + 1
          e%ounit = e%munit
        end if
      end if
      return ! Done
    end subroutine done

    subroutine col_head
      if (e%lentry == -3) then ! Take care of previous
        k = lbuf - lbufh
        if (k > hlen ) then
          call mess_error(1, "$H text has insufficient room for heading")
          return
        else
          if (k == 0) return
        end if
        select case(hj)
        case ('C') ! Center
          j = (hlen - k + 1) / 2
          buf(lenbuf-k+1:lenbuf) = buf(lbufh+1:lbuf) ! A copy
          buf(lbufh+1:lbufh+j) = ' ' ! Blanks at the start
          buf(lbufh+j+1:lbufh+j+k) = buf(lenbuf-k+1:lenbuf)
          buf(lbufh+j+k+1:lbuf+hlen) = ' ' ! Blanks at the end
        case ('R') ! Right justify
          buf(lenbuf-k+1:lenbuf) = buf(lbufh+1:lbuf) ! A copy
          buf(lbufh+1:lbufh+hlen-k) = ' ' ! Blanks at the start
          buf(lbufh+hlen -k+1:lbufh+hlen) = buf(lenbuf-k+1:lenbuf)
        case ('L') ! Left justify
          buf(lbuf+1:lbufh+hlen) = ' '
        case default
          c = hj
          call mess_error(11, "C, R, or L (upper case) must be used&
            & for column justification and text terminated with&
            & the special character.  You have a ")
          return
        end select
        lbuf = lbufh + hlen
        lbufh = lbuf
        if (itb == 1) then
          e%tbrk(1) = lbuf
          itb = 2
        end if
      end if
      if (lbuf <= e%lenlin) then
        e%tbrk(itb) = lbuf
        e%brk_lines = itb - 2
      end if
      e%lentry = -3
      lentxt = ntext + index(text(ntext+1:),csh)
      if (lentxt <= ntext+2) then ! No room for column head
        if (lentxt == lentxth) then ! Reached end of the text
          do while (lbuf > e%lenlin) ! Flush what we can
            call bufout ! Not >= above so next e%tbrk is set.
          end do
          lbufh = lbuf
          e%tbrk(itb) = lbuf
          e%brk_lines = itb - 2
          ntext = lentxt
          return
        else
          lentxt = lentxth
          call mess_error(1, &
            & "$H text does not end with character following the $H, or&
            & no space for heading info. after the special character.")
        end if
        return
      end if
      hlen = getint()
      if (hlen == 0) then
        hlen = e%kdf + 8
        if (text(ntext+1:ntext+1) == '0') then
          hlen = 2*hlen + 3 ! Assume complex  " (num,num)"
          ntext = ntext + 1
        end if
      end if
      ntext = ntext + 1
      hj = text(ntext:ntext)
    end subroutine col_head


    subroutine finish_table
      do while (lbuf > 0)
        call  bufout
      end do
      if (e%brk_lines > 0) then
        endfile e%sunit ! The Cray ftn compiler appears to ignore this statment.
        do i = 1, e%brk_lines
          rewind(e%sunit)
          do j = 1, i - 1 ! END below should not happen
            read (e%sunit, '(A)', END=100) buf
          end do
          do
            read (e%sunit, '(A)', END = 100) buf(1:e%tbrk(i+2))
            write (e%munit, '(A)') buf(1:e%tbrk(i+2))
            do j = 1, e%brk_lines - 1
              read (e%sunit, '(A)', END = 100) buf
            end do
          end do
100       continue
        end do
        e%brk_lines = 0
        close(e%sunit)
        e%sunit = 0
      end if
      return
    end subroutine finish_table

    subroutine mess_error(ierr, errtxt)
      character (len=*), intent (in) :: errtxt
      integer, intent(in) :: ierr
      integer :: i, j
      i = 1
      j = min(ntext+2, len(text))
      if (j - i > e%lenlin) i = j - e%lenlin
      k = e%lentry
      e%lentry = 0
      if (lbuf /= 0) write(e%eunit, '(//"Was getting ready to output:"/A)')&
        buf(1:lbuf)
      call messy(e, etext)
      write (e%eunit, '(''"'',A,''"'')') text(i:j)
      write (e%eunit, '("All of text"/A)') text
      select case (ierr)
      case (1)
        write (e%eunit, '(A)') errtxt
      case (2)
        write (e%eunit, '("$",A," used, but the call does not contain ",A)')&
          & errtxt(1:1), errtxt(2:)
      case (3)

      case (4)
        write (e%eunit, '("Attempt to access text past the end")')
      case (5)
        write (e%eunit,&
          &'("Requested data is above upper bound ", I0," of ",A)')&
          & errdat, errtxt
      case (6)
        write (e%eunit, '("$Z",A," requests data beyond the size of ",A)')&
          & errtxt(1:1), errtxt(2:)
      case (7)
        write (e%eunit, '(A," used with no format defined")') errtxt
      case (8)
        write (e%eunit, '("Sizes of idat and rdat are ",i0, " and ", i0,&
        & "."/"They must be equal for sparse output.")') bndi, bndr
      case (9)
        write (e%eunit, '("Bit strings require ", i0, " words, but the space&
          & n idat for this vector of bit strings" /&
          & " is", i0, " which is not a multiple of ", i0, ".")') nwords,&
          & bndi - nidat, nwords
      case (10)
        write(e%eunit, '("Current numbers require a line length of ", i0,&
        & " but your line length is only ", i0,"."/&
        & "Increase the line length or decrease the number of digits.")')&
        & neg, e%lenlin
      case (11)
        write (e%eunit, '(A, A1,".")') errtxt, c
      end select
      call messy(e, "$B")
      e%lentry = k
      e%ounit = e%munit
      e%errcnt = e%errcnt + 1000000
      ntext = lentxt + 1 ! So we return from outer calls
      return
    end subroutine mess_error
  end subroutine messy

