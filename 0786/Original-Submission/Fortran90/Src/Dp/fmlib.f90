
!     FM 1.1                  David M. Smith               5-19-97



!  The FM routines in this package perform floating-point
!  multiple-precision arithmetic, and the IM routines perform
!  integer multiple-precision arithmetic.



!  1. INITIALIZING THE PACKAGE

!  Before calling any routine in the package, several variables in
!  the common blocks /FMUSER/, /FM/, /FMBUFF/, and /FMSAVE/ must be
!  initialized.  These four common blocks contain information that
!  is saved between calls, so they should be declared in the main
!  program.

!  Subroutine FMSET initializes these variables to default values and
!  defines all machine-dependent values in the package.  After calling
!  FMSET once at the start of a program, the user may sometimes want
!  to reset some of the variables in these common blocks.  These
!  variables are described below.



!  2.  REPRESENTATION OF FM NUMBERS

!  MBASE is the base in which the arithmetic is done.  MBASE must be
!        bigger than one, and less than or equal to the square root of
!        the largest representable integer.  For best efficiency MBASE
!        should be large, but no more than about 1/4 of the square root
!        of the largest representable integer.  Input and output
!        conversions are much faster when MBASE is a power of ten.

!  NDIG  is the number of base MBASE digits that are carried in the
!        multiple precision numbers.  NDIG must be at least two.  The
!        upper limit for NDIG is defined in the PARAMETER statement at
!        the top of each routine and is restricted only by the amount
!        of memory available.

!  Sometimes it is useful to dynamically vary NDIG during the program.
!  Use FMEQU to round numbers to lower precision or zero-pad them to
!  higher precision when changing NDIG.

!  It is rare to need to change MBASE during a program.  Use FMCONS to
!  reset some saved constants that depend on MBASE.  FMCONS should be
!  called immediately after changing MBASE.

!  There are two representations for a floating multiple precision
!  number.  The unpacked representation used by the routines while
!  doing the computations is base MBASE and is stored in NDIG+2 words.
!  A packed representation is available to store the numbers in the
!  user's program in compressed form.  In this format, the NDIG
!  (base MBASE) digits of the mantissa are packed two per word to
!  conserve storage.  Thus the external, packed form of a number
!  requires (NDIG+1)/2+2 words.

!  This version uses double precision arrays to hold the numbers.
!  Version 1.0 of FM used integer arrays, which are faster on some
!  machines.  The package can easily be changed to use integer
!  arrays -- see section 11 on EFFICIENCY below.

!  The unpacked format of a floating multiple precision number is as
!  follows.  A number MA is kept in an array with MA(1) containing
!  the exponent and MA(2) through MA(NDIG+1) containing one digit of
!  the mantissa, expressed in base MBASE.  The array is dimensioned
!  to start at MA(0), with the approximate number of bits of precision
!  stored in MA(0).  This precision value is intended to be used by FM
!  functions that need to monitor cancellation error in addition and
!  subtraction. The cancellation monitor code is usually disabled for
!  user calls, and FM functions only check for cancellation when they
!  must.  Tracking cancellation causes most routines to run slower,
!  with addition and subtraction being affected the most.

!  The exponent is a power of MBASE and the implied radix point is
!  immediately before the first digit of the mantissa.  Every nonzero
!  number is normalized so that the second array element (the first
!  digit of the mantissa) is nonzero.

!  In both representations the sign of the number is carried on the
!  second array element only.  Elements 3,4,... are always nonnegative.
!  The exponent is a signed integer and may be as large in magnitude as
!  MXEXP (defined in FMSET).

!  For MBASE = 10,000 and NDIG = 4, the number -pi would have these
!  representations:
!                   Word 1         2         3         4         5

!         Unpacked:      1        -3      1415      9265      3590
!         Packed:        1    -31415  92653590

!  Word 0 would be 42 in both formats, indicating that the mantissa
!  has about 42 bits of precision.

!  Because of normalization in a large base, the equivalent number
!  of base 10 significant digits for an FM number may be as small as
!  LOG10(MBASE)*(NDIG-1) + 1.

!  The integer routines use the FMLIB format to represent numbers,
!  without the number of digits (NDIG) being fixed.  Integers in IM
!  format are essentially variable precision, using the minimum number
!  of words to represent each value.

!  For programs using both FM and IM numbers, FM routines should not
!  be called with IM numbers, and IM routines should not be called
!  with FM numbers, since the implied value of NDIG used for an IM
!  number may not match the explicit NDIG expected by an FM routine.
!  Use the conversion routines IMFM2I and IMI2FM to change between
!  the FM and IM formats.



!  3. INPUT/OUTPUT ROUTINES

!  All versions of the input routines perform free-format conversion
!  from characters to FM numbers.

!  a. Conversion to or from a character array

!     FMINP converts from a character*1 array to an FM number.

!     FMOUT converts an FM number to base 10 and formats it for output
!           as an array of type character*1.  The output is left
!           justified in the array, and the format is defined by two
!           variables in common, so that a separate format definition
!           does not have to be provided for each output call.

!     The user sets JFORM1 and JFORM2 to determine the output format.

!     JFORM1 = 0     E   format       ( .314159M+6 )
!            = 1     1PE format       ( 3.14159M+5 )
!            = 2     F   format       ( 314159.000 )

!     JFORM2 is the number of significant digits to display (if
!            JFORM1 = 0 or 1).  If JFORM2.EQ.0 then a default number
!            of digits is chosen.  The default is roughly the full
!            precision of the number.
!     JFORM2 is the number of digits after the decimal point (if
!            JFORM1 = 2).  See the FMOUT documentation for more details.

!  b. Conversion to or from a character string

!     FMST2M converts from a character string to an FM number.

!     FMFORM converts an FM number to a character string according to
!            a format provided in each call.  The format description
!            is more like that of a Fortran FORMAT statement, and
!            integer or fixed-point output is right justified.

!  c. Direct read or write

!     FMPRNT uses FMOUT to print one FM number.

!     FMFPRT uses FMFORM to print one FM number.

!     FMWRIT writes FM numbers for later input using FMREAD.

!     FMREAD reads FM numbers written by FMWRIT.

!  The values given to JFORM1 and JFORM2 can be used to define a
!  default output format when FMOUT or FMPRNT are called.  The
!  explicit format used in a call to FMFORM or FMFPRT overrides
!  the settings of JFORM1 and JFORM2.

!  KW is the unit number to be used for standard output from
!     the package, including error and warning messages, and
!     trace output.

!  For multiple precision integers, the corresponding routines
!  IMINP, IMOUT, IMST2M, IMFORM, IMPRNT, IMFPRT, IMWRIT, and
!  IMREAD provide similar input and output conversions.  For
!  output of IM numbers, JFORM1 and JFORM2 are ignored and
!  integer format (JFORM1=2, JFORM2=0) is used.

!  For further description of these routines, see sections
!  9 and 10 below.



!  4. ARITHMETIC TRACING

!  NTRACE and LVLTRC control trace printout from the package.

!  NTRACE =  0   No printout except warnings and errors.
!         =  1   The result of each call to one of the routines
!                   is printed in base 10, using FMOUT.
!         = -1   The result of each call to one of the routines
!                   is printed in internal base MBASE format.
!         =  2   The input arguments and result of each call to one
!                   of the routines is printed in base 10, using FMOUT.
!         = -2   The input arguments and result of each call to one
!                   of the routines is printed in base MBASE format.

!  LVLTRC defines the call level to which the trace is done.  LVLTRC = 1
!         means only FM routines called directly by the user are traced,
!         LVLTRC = 2 also prints traces for FM routines called by other
!         FM routines called directly by the user, etc.

!  In the above description, internal MBASE format means the number is
!  printed as it appears in the array --- an exponent followed by NDIG
!  base MBASE digits.



!  5. ERROR CONDITIONS

!  KFLAG is a condition parameter returned by the package after each
!        call to one of the routines.  Negative values indicate
!        conditions for which a warning message will be printed
!        unless KWARN = 0.  Positive values indicate conditions
!        that may be of interest but are not errors.
!        No warning message is printed if KFLAG is nonnegative.

!    KFLAG =  0     Normal operation.

!          =  1     One of the operands in FMADD or FMSUB was
!                       insignificant with respect to the other, so
!                       that the result was equal to the argument of
!                       larger magnitude.
!          =  2     In converting an FM number to a one word integer
!                       in FMM2I, the FM number was not exactly an
!                       integer.  The next integer toward zero was
!                       returned.

!          = -1     NDIG was less than 2 or more than NDIGMX.
!          = -2     MBASE was less than 2 or more than MXBASE.
!          = -3     An exponent was out of range.
!          = -4     Invalid input argument(s) to an FM routine.
!                        UNKNOWN was returned.
!          = -5     + or - OVERFLOW was generated as a result from an
!                        FM routine.
!          = -6     + or - UNDERFLOW was generated as a result from an
!                        FM routine.
!          = -7     The input string (array) to FMINP was not legal.
!          = -8     The character array was not large enough in an
!                   input or output routine.
!          = -9     Precision could not be raised enough to provide all
!                        requested guard digits.  Increasing NDIGMX in
!                        all the PARAMETER statements may fix this.
!                        UNKNOWN was returned.
!          = -10    An FM input argument was too small in magnitude to
!                        convert to the machine's single or double
!                        precision in FMM2SP or FMM2DP.  Check that the
!                        definitions of SPMAX and DPMAX in FMSET are
!                        correct for the current machine.
!                        Zero was returned.

!  When a negative KFLAG condition is encountered, the value of KWARN
!  determines the action to be taken.

!  KWARN = 0     Execution continues and no message is printed.
!        = 1     A warning message is printed and execution continues.
!        = 2     A warning message is printed and execution stops.

!  The default setting is KWARN = 1.

!  When an overflow or underflow is generated for an operation in which
!  an input argument was already an overflow or underflow, no additional
!  message is printed.  When an unknown result is generated and an input
!  argument was already unknown, no additional message is printed.  In
!  these cases the negative KFLAG value is still returned.

!  IM routines handle exceptions like OVERFLOW or UNKNOWN in the same
!  way as FM routines.  When using IMMPY, the product of two large
!  positive integers will return +OVERFLOW.  The routine IMMPYM can
!  be used to obtain a modular result without overflow.  The largest
!  representable IM integer is MBASE**NDIGMX - 1.  For example, if
!  MBASE is 10**7 and NDIGMX is set to 256, integers less than 10**1792
!  can be used.



!  6. OTHER PARAMETERS

!  KRAD = 0     All angles in the trigonometric functions and
!                  inverse functions are measured in degrees.
!       = 1     All angles are measured in radians. (Default)


!  KROUND = 0   All final results are chopped (rounded toward
!                  zero).  Intermediate results are rounded.
!         = 1   All results are rounded to the nearest FM
!                  number, or to the value with an even last
!                  digit if the result is halfway between two
!                  FM numbers. (Default)

!  KSWIDE defines the maximum screen width to be used for
!         all unit KW output.  Default is 80.

!  KESWCH controls the action taken in FMINP and other input routines
!         for strings like 'E7' that have no digits before the exponent
!         field.  Default is for 'E7' to translate like '1.0E+7'.

!  CMCHAR defines the exponent letter to be used for FM variable
!         output.  Default is 'M', as in 1.2345M+678.

!  KDEBUG = 0   Error checking is not done for valid input arguments
!               and parameters like NDIG and MBASE upon entry to
!               each routine. (Default)
!         = 1   Some error checking is done.  (Slower speed)

!  See FMSET for additional description of these and other variables
!  defining various FM conditions.



!  7. ARRAY DIMENSIONS

!  The dimensions of the arrays in the FM package are defined using
!  a PARAMETER statement at the top of each routine.  The size of
!  these arrays depends on the values of parameters NDIGMX and NBITS.
!  NDIGMX is the maximum value the user may set for NDIG.
!  NBITS  is the number of bits used to represent integers for a
!         given machine.  See the EFFICIENCY discussion below.

!  The standard version of FMLIB sets NDIGMX = 256, so on a 32-bit
!  machine using MBASE = 10**7 the maximum precision is about
!  7*255+1 = 1786 significant digits.  To change dimensions so that
!  10,000 significant digit calculation can be done, NDIGMX needs to
!  be at least  10**4/7 + 5 = 1434.  This allows for a few user guard
!  digits to be defined when the package is initialized using
!  CALL FMSET(10000).  Changing 'NDIGMX = 256' to 'NDIGMX = 1434'
!  everywhere in the package and the user's calling program will
!  define all the new array sizes.

!  If NDIG much greater than 256 is to be used and elementary functions
!  will be needed, they will be faster if array MJSUMS is larger. The
!  parameter defining the size of MJSUMS is set in the standard version
!  by LJSUMS = 8*(LUNPCK+2).  The 8 means that up to eight concurrent
!  sums can be used by the elementary functions.  The approximate number
!  needed for best speed is given by the formula
!      0.051*Log(MBASE)*NDIG**(1/3) + 1.85
!  For example, with MBASE=10**7 and NDIG=1434 this gives 11.  Changing
!  'LJSUMS = 8*(LUNPCK+2)' to 'LJSUMS =11*(LUNPCK+2)' everywhere in the
!  package and the user's calling program will give slightly better
!  speed.

!  FM numbers in packed format have dimension 0:LPACK, and those
!  in unpacked format have dimension 0:LUNPCK.



!  8. PORTABILITY

!  In FMSET there is some machine-dependent code that attempts to
!  approximate the largest representable integer value.  The current
!  code works on all machines tested, but if an FM run fails, check
!  the MAXINT and INTMAX loops in FMSET.  Values for SPMAX and DPMAX
!  are also defined in FMSET that should be set to values near overflow
!  for single precision and double precision.  Setting KDEBUG = 1 may
!  also identify some errors if a run fails.

!  Some compilers object to a function like FMCOMP with side effects
!  such as changing KFLAG or other common variables.  Blocks of code
!  in FMCOMP and IMCOMP that modify common are identified so they may
!  be removed or commented out to produce a function without side
!  effects.  This disables trace printing in FMCOMP and IMCOMP, and
!  error codes are not returned in KFLAG.  See FMCOMP and IMCOMP for
!  further details.


!  9.  LIST OF ROUTINES

!  These are the FM routines that are designed to be called by
!  the user.  All are subroutines except logical function FMCOMP.
!  MA, MB, MC refer to FM format numbers.

!  In each case it is permissible to use the same array more than
!  once in the calling sequence.  The statement MA = MA*MA can
!  be written CALL FMMPY(MA,MA,MA).

!  For each of these routines there is also a version available for
!  which the argument list is the same but all FM numbers are in packed
!  format.  The routines using packed numbers have the same names except
!  'FM' is replaced by 'FP' at the start of each name.


!  FMABS(MA,MB)         MB = ABS(MA)

!  FMACOS(MA,MB)        MB = ACOS(MA)

!  FMADD(MA,MB,MC)      MC = MA + MB

!  FMADDI(MA,IVAL)      MA = MA + IVAL   Increment an FM number by a one
!                                        word integer.  Note this call
!                                        does not have an "MB" result
!                                        like FMDIVI and FMMPYI.

!  FMASIN(MA,MB)        MB = ASIN(MA)

!  FMATAN(MA,MB)        MB = ATAN(MA)

!  FMATN2(MA,MB,MC)     MC = ATAN2(MA,MB)

!  FMBIG(MA)            MA = Biggest FM number less than overflow.

!  FMCHSH(MA,MB,MC)     MB = COSH(MA),  MC = SINH(MA).  Faster than
!                            making two separate calls.

!  FMCOMP(MA,LREL,MB)        Logical comparison of MA and MB.
!                            LREL is a CHARACTER*2 value identifying
!                            which comparison is made.
!                            Example:  IF (FMCOMP(MA,'GE',MB)) ...

!  FMCONS                    Set several saved constants that depend
!                            on MBASE, the base being used.  FMCONS
!                            should be called immediately after
!                            changing MBASE.

!  FMCOS(MA,MB)         MB = COS(MA)

!  FMCOSH(MA,MB)        MB = COSH(MA)

!  FMCSSN(MA,MB,MC)     MB = COS(MA),  MC = SIN(MA).  Faster than
!                            making two separate calls.

!  FMDIG(NSTACK,KST)         Find a set of precisions to use during
!                            Newton iteration for finding a simple
!                            root starting with about double
!                            precision accuracy.

!  FMDIM(MA,MB,MC)      MC = DIM(MA,MB)

!  FMDIV(MA,MB,MC)      MC = MA/MB

!  FMDIVI(MA,IVAL,MB)   MB = MA/IVAL   IVAL is a one word integer.

!  FMDP2M(X,MA)         MA = X    Convert from double precision to FM.

!  FMDPM(X,MA)          MA = X    Convert from double precision to FM.
!                                 Much faster than FMDP2M, but MA agrees
!                                 with X only to D.P. accuracy.  See
!                                 the comments in the two routines.

!  FMEQ(MA,MB)          MB = MA   Both have precision NDIG.
!                                 This is the version to use for
!                                 standard  B = A  statements.

!  FMEQU(MA,MB,NA,NB)   MB = MA   Version for changing precision.
!                                 MA has NA digits (i.e., MA was
!                                 computed using NDIG = NA), and MB
!                                 will be defined having NB digits.
!                                 MB is zero-padded if NB.GT.NA
!                                 MB is rounded if NB.LT.NA

!  FMEXP(MA,MB)         MB = EXP(MA)

!  FMFORM(FORM,MA,STRING)    MA is converted to a character string
!                               using format FORM and returned in
!                               STRING.  FORM can represent I, F,
!                               E, or 1PE formats.  Example:
!                               CALL FMFORM('F60.40',MA,STRING)

!  FMFPRT(FORM,MA)           Print MA on unit KW using FORM format.

!  FMI2M(IVAL,MA)       MA = IVAL   Convert from one word integer
!                                   to FM.

!  FMINP(LINE,MA,LA,LB) MA = LINE   Input conversion.
!                                   Convert LINE(LA) through LINE(LB)
!                                   from characters to FM.

!  FMINT(MA,MB)         MB = INT(MA)    Integer part of MA.

!  FMIPWR(MA,IVAL,MB)   MB = MA**IVAL   Raise an FM number to a one
!                                       word integer power.

!  FMLG10(MA,MB)        MB = LOG10(MA)

!  FMLN(MA,MB)          MB = LOG(MA)

!  FMLNI(IVAL,MA)       MA = LOG(IVAL)   Natural log of a one word
!                                        integer.

!  FMM2DP(MA,X)         X  = MA     Convert from FM to double precision.

!  FMM2I(MA,IVAL)       IVAL = MA   Convert from FM to integer.

!  FMM2SP(MA,X)         X  = MA     Convert from FM to single precision.

!  FMMAX(MA,MB,MC)      MC = MAX(MA,MB)

!  FMMIN(MA,MB,MC)      MC = MIN(MA,MB)

!  FMMOD(MA,MB,MC)      MC = MA mod MB

!  FMMPY(MA,MB,MC)      MC = MA*MB

!  FMMPYI(MA,IVAL,MB)   MB = MA*IVAL    Multiply by a one word integer.

!  FMNINT(MA,MB)        MB = NINT(MA)   Nearest FM integer.

!  FMOUT(MA,LINE,LB)    LINE = MA   Convert from FM to character.
!                                   LINE is a character array of
!                                   length LB.

!  FMPI(MA)             MA = pi

!  FMPRNT(MA)                Print MA on unit KW using current format.

!  FMPWR(MA,MB,MC)      MC = MA**MB

!  FMREAD(KREAD,MA)     MA   is returned after reading one (possibly
!                            multi-line) FM number on unit KREAD.  This
!                            routine reads numbers written by FMWRIT.

!  FMRPWR(MA,K,J,MB)    MB = MA**(K/J)  Rational power.  Faster than
!                            FMPWR for functions like the cube root.

!  FMSET(NPREC)              Set default values and machine-dependent
!                            variables to give at least NPREC base 10
!                            digits plus three base 10 guard digits.
!                            Must be called to initialize FM package.

!  FMSIGN(MA,MB,MC)     MC = SIGN(MA,MB)   Sign transfer.

!  FMSIN(MA,MB)         MB = SIN(MA)

!  FMSINH(MA,MB)        MB = SINH(MA)

!  FMSP2M(X,MA)         MA = X   Convert from single precision to FM.

!  FMSQR(MA,MB)         MB = MA*MA   Faster than FMMPY.

!  FMSQRT(MA,MB)        MB = SQRT(MA)

!  FMST2M(STRING,MA)    MA = STRING
!                            Convert from character string to FM.
!                            Often more convenient than FMINP, which
!                            converts an array of CHARACTER*1 values.
!                            Example:   CALL FMST2M('123.4',MA).

!  FMSUB(MA,MB,MC)      MC = MA - MB

!  FMTAN(MA,MB)         MB = TAN(MA)

!  FMTANH(MA,MB)        MB = TANH(MA)

!  FMULP(MA,MB)         MB = One Unit in the Last Place of MA.

!  FMWRIT(KWRITE,MA)         Write MA on unit KWRITE.
!                            Multi-line numbers will have '&' as the
!                            last nonblank character on all but the last
!                            line.  These numbers can then be read
!                            easily using FMREAD.


!  These are the integer routines that are designed to be called by
!  the user.  All are subroutines except logical function IMCOMP.
!  MA, MB, MC refer to IM format numbers.  In each case the version
!  of the routine to handle packed IM numbers has the same name,
!  with 'IM' replaced by 'IP'.

!  IMABS(MA,MB)         MB = ABS(MA)

!  IMADD(MA,MB,MC)      MC = MA + MB

!  IMBIG(MA)            MA = Biggest IM number less than overflow.

!  IMCOMP(MA,LREL,MB)        Logical comparison of MA and MB.
!                            LREL is a CHARACTER*2 value identifying
!                            which comparison is made.
!                            Example:  IF (IMCOMP(MA,'GE',MB)) ...

!  IMDIM(MA,MB,MC)      MC = DIM(MA,MB)

!  IMDIV(MA,MB,MC)      MC = int(MA/MB)
!                            Use IMDIVR if the remainder is also needed.

!  IMDIVI(MA,IVAL,MB)   MB = int(MA/IVAL)
!                            IVAL is a one word integer.  Use IMDVIR
!                            to get the remainder also.

!  IMDIVR(MA,MB,MC,MD)  MC = int(MA/MB),   MD = MA mod MB
!                            When both the quotient and remainder are
!                            needed, this routine is twice as fast as
!                            calling both IMDIV and IMMOD.

!  IMDVIR(MA,IVAL,MB,IREM)   MB = int(MA/IVAL),   IREM = MA mod IVAL
!                            IVAL and IREM are one word integers.

!  IMEQ(MA,MB)          MB = MA

!  IMFM2I(MAFM,MB)      MB = MAFM  Convert from real (FM) format
!                                  to integer (IM) format.

!  IMFORM(FORM,MA,STRING)    MA is converted to a character string
!                               using format FORM and returned in
!                               STRING.  FORM can represent I, F,
!                               E, or 1PE formats.  Example:
!                               CALL IMFORM('I70',MA,STRING)

!  IMFPRT(FORM,MA)           Print MA on unit KW using FORM format.

!  IMGCD(MA,MB,MC)      MC = greatest common divisor of MA and MB.

!  IMI2FM(MA,MBFM)    MBFM = MA  Convert from integer (IM) format
!                                to real (FM) format.

!  IMI2M(IVAL,MA)       MA = IVAL   Convert from one word integer
!                                   to IM.

!  IMINP(LINE,MA,LA,LB) MA = LINE   Input conversion.
!                                   Convert LINE(LA) through LINE(LB)
!                                   from characters to IM.

!  IMM2DP(MA,X)         X  = MA     Convert from IM to double precision.

!  IMM2I(MA,IVAL)       IVAL = MA   Convert from IM to one word integer.

!  IMMAX(MA,MB,MC)      MC = MAX(MA,MB)

!  IMMIN(MA,MB,MC)      MC = MIN(MA,MB)

!  IMMOD(MA,MB,MC)      MC = MA mod MB

!  IMMPY(MA,MB,MC)      MC = MA*MB

!  IMMPYI(MA,IVAL,MB)   MB = MA*IVAL    Multiply by a one word integer.

!  IMMPYM(MA,MB,MC,MD)  MD = MA*MB mod MC
!                            Slightly faster than calling IMMPY and
!                            IMMOD separately, and it works for cases
!                            where IMMPY would return OVERFLOW.

!  IMOUT(MA,LINE,LB)    LINE = MA   Convert from IM to character.
!                                   LINE is a character array of
!                                   length LB.

!  IMPMOD(MA,MB,MC,MD)       MD = MA**MB mod MC

!  IMPRNT(MA)                Print MA on unit KW.

!  IMPWR(MA,MB,MC)      MC = MA**MB

!  IMREAD(KREAD,MA)     MA   is returned after reading one (possibly
!                            multi-line) IM number on unit KREAD.  This
!                            routine reads numbers written by IMWRIT.

!  IMSIGN(MA,MB,MC)     MC = SIGN(MA,MB)   Sign transfer.

!  IMSQR(MA,MB)         MB = MA*MA   Faster than IMMPY.

!  IMST2M(STRING,MA)    MA = STRING
!                            Convert from character string to IM.
!                            Often more convenient than IMINP, which
!                            converts an array of CHARACTER*1 values.
!                            Example:   CALL IMST2M('12345678901',MA).

!  IMSUB(MA,MB,MC)      MC = MA - MB

!  IMWRIT(KWRITE,MA)         Write MA on unit KWRITE.
!                            Multi-line numbers will have '&' as the
!                            last nonblank character on all but the last
!                            line.  These numbers can then be read
!                            easily using IMREAD.

!  Many of the IM routines call FM routines, but none of the FM
!  routines call IM routines, so the IM routines can be omitted
!  if none are called explicitly from a program.



!  10. NEW FOR VERSION 1.1

!  Version 1.0 used integer arrays and integer arithmetic internally
!  to perform the multiple precision operations.  Version 1.1 uses
!  double precision arithmetic and arrays internally.  This is usually
!  faster at higher precisions, and on many machines it is also faster
!  at lower precisions.  Version 1.1 is written so that the arithmetic
!  used can easily be changed from double precision to integer, or any
!  other available arithmetic type.  This permits the user to make the
!  best use of a given machine's arithmetic hardware.
!  See the EFFICIENCY discussion below.

!  Several routines have undergone minor modification, but only a few
!  changes should affect programs that used FM 1.0.  Many of the
!  routines are faster in version 1.1, because code has been added to
!  take advantage of special cases for individual functions instead of
!  using general formulas that are more compact.  For example, there
!  are separate routines using series for SINH and COSH instead of
!  just calling EXP.

!  FMEQU was the only routine that required the user to give the value
!        of the current precision.  This was to allow automatic
!        rounding or zero-padding when changing precision.  Since few
!        user calls change precision, a new routine has been added for
!        this case.
!        FMEQ now handles this case and has a simple argument list that
!             does not include the value of NDIG.
!        FMEQU is used for changing precision.
!        See the list of FM routines above for details.

!  All variable names beginning with M in the package are now declared
!  as double precision, so FM common blocks in the user's program need
!  D.P. declarations, and FM variables (arrays) used in the calling
!  program need to be D.P.

!  /FMUSER/ is a common block holding parameters that define the
!           arithmetic to be used and other user options.  Several
!           new variables have been added, including screen width to
!           be used for output.  See above for further description.

!  /FMSAVE/ is a common block for saving constants to avoid
!           re-computing them.  Several new variables have been added.

!  /FMBUFF/ is a common block containing a character array used to
!           format FM numbers for output.  Two new items have been
!           added.

!  New routines:

!  All the IM routines are new for version 1.1.

!  FMADDI increments an FM number by a small integer.
!         It runs in O(1) time, on the average.

!  FMCHSH returns both SINH(MA) and COSH(MA).
!         When both are needed, this is almost twice as fast
!         as making separate calls to FMCOSH and FMSINH.

!  FMCSSN returns both SIN(MA) and COS(MA).
!         When both are needed, this is almost twice as fast
!         as making separate calls to FMCOS and FMSIN.

!  FMFORM uses a format string to convert an FM number to a
!         character string.

!  FMFPRT prints an FM number using a format string.

!  FMREAD reads an FM number written using FMWRIT.

!  FMRPWR computes an FM number raised to a rational power.  For cube
!         roots and similar rational powers it is usually much faster
!         than FMPWR.

!  FMSQR  squares an FM number.  It is faster than using FMMPY.

!  FMST2M converts character strings to FM format.  Since FMINP converts
!         character arrays, this routine can be more convenient for
!         easily defining an FM number.
!         For example, CALL FMST2M('123.4',MA).

!  FMWRIT writes an FM number using a format for multi-line numbers
!         with '&' at the end of all but the last line of a multi-line
!         number.  This allows automatic reading of FM numbers without
!         needing to know the base, precision or format under which they
!         were written.

!  One extra word has been added to the dimensions of all FM numbers.
!  Word zero in each array contains a value used to monitor cancellation
!  error arising from addition or subtraction.  This value approximates
!  the number of bits of precision for an FM value.  It allows higher
!  level FM functions to detect cases where too much cancellation has
!  occurred.  KACCSW is a switch variable in COMMON /FM/ used internally
!  to enable cancellation error monitoring.


!  11. EFFICIENCY

!  To take advantage of hardware architecture on different machines, the
!  package has been designed so that the arithmetic used to perform the
!  multiple precision operations can easily be changed.  All variables
!  that must be changed to get a different arithmetic have names
!  beginning with 'M' and are declared using REAL (KIND(0.0D0)) :: m....

!  For example, to change the package to use integer arithmetic
!  internally, make these two changes everywhere in the package:
!  change  'REAL (KIND(0.0D0)) :: m'  to  'INTEGER m',
!  change  'DINT('  to  'INT('.
!  On some systems, changing  'DINT('  to  '('  may give better speed.

!  When changing to a different type of arithmetic, all FM common blocks
!  and arrays in the user's program must be changed to agree.  In a few
!  places in FM, where a DINT function is not supposed to be changed, it
!  is spelled 'DINT (' so the global change will not find it.

!  This version restricts the base used to be also representable in
!  integer variables, so using precision above double usually does not
!  save much time unless integers can also be declared at a higher
!  precision.  Using IEEE Extended would allow a base of around 10**9
!  to be chosen, but the delayed digit-normalization method used for
!  multiplication and division means that a slightly smaller base like
!  10**8 would usually run faster.  This would usually not be much
!  faster than using 10**7 with double precision.

!  The value of NBITS defined as a parameter in most FM routines
!  refers to the number of bits used to represent integers in an
!  M-variable word.  Typical values for NBITS are:  24 for IEEE single
!  precision, 32 for integer, 53 for IEEE double precision.  NBITS
!  controls only array size, so setting it too high is ok, but then
!  the program will use more memory than necessary.

!  For cases where special compiler directives or minor re-writing
!  of the code may improve speed, several of the most important
!  loops in FM are identified by comments containing the string
!  '(Inner Loop)'.

! --------------------------------------------------------------------
! --------------------------------------------------------------------


    SUBROUTINE fmset(nprec)

!  Initialize the values in common that must be set before calling
!  other FM routines.

!  Base and precision will be set to give at least NPREC+3 decimal
!  digits of precision (giving the user three base ten guard digits).

!  MBASE is set to a large power of ten.
!  JFORM1 and JFORM2 are set to 1PE format displaying NPREC
!  significant digits.

!  The trace option is set off.
!  The mode for angles in trig functions is set to radians.
!  The rounding mode is set to symmetric rounding.
!  Warning error message level is set to 1.
!  Cancellation error monitor is set off.
!  Screen width for output is set to 80 columns.
!  The exponent character for FM output is set to 'M'.
!  Debug error checking is set off.

!  KW, the unit number for all FM output, is set to 6.

!  The size of all arrays is controlled by defining two parameters:
!  NDIGMX is the maximum value the user can set NDIG,
!  NBITS  is the number of bits used to represent integers in an
!         M-variable word.

      IMPLICIT NONE

!             Define the array sizes:

!  Here are all the common blocks used in FM.

!  /FMUSER/, /FM/, /FMBUFF/, and /FMSAVE/ should also be declared in the
!  main program, because some compilers allocate and free space used for
!  labelled common that is declared only in subprograms.  This causes
!  the saved information to be lost.

!             FMUSER contains values that may need to be
!                    changed by the calling program.

!             FM contains the work array used by the low-level
!                arithmetic routines, definitions for overflow
!                and underflow thresholds, and other
!                machine-dependent values.

!             FMSAVE contains information about saved constants.

!             MJSUMS is an array that can contain several FM numbers
!             being used to accumulate concurrent sums in exponential
!             and trigonometric functions.  When NDIGMX = 256, eight is
!             about the maximum number of sums needed (but this depends
!             on MBASE).  For larger NDIGMX, dimensioning MJSUMS to hold
!             more than eight FM numbers could increase the speed of the
!             functions.

!             FMWA contains two work arrays similar to MWA.  They are
!             used in routines FMDIVD, FMMPYD, and FMMPYE.

!             CMBUFF is a character array used by FMPRNT for printing
!                    output from FMOUT.  This array may also be used
!                    for calls to FMOUT from outside the FM package.
!             CMCHAR is the letter used before the exponent field
!                    in FMOUT.  It is defined in FMSET.
!             NAMEST is a stack for names of the routines.  It is
!                    used for trace printing and error messages.

!             FM1 contains scratch arrays for temporary storage of FM
!             numbers while computing various functions.

!             FMPCK contains scratch arrays used to hold input arguments
!             in unpacked format when the packed versions of functions
!             are used.

! .. Intrinsic Functions ..
      INTRINSIC dble, ichar, int, log, log10, max, min, sqrt
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: nprec
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ml, mld2, mlm1
      REAL (KIND(0.0D0)) :: one, temp, two, yt
      INTEGER :: j, k, kpt, l, npsave
! ..
! .. Local Arrays ..
      INTEGER :: ltypes(21), lvals(21)
      CHARACTER (1) :: lchars(21)
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmdbl, fmmset
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mjsums(0:ljsums), &
        mlbsav(0:lunpck), mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), &
        mln4(0:lunpck), mpa(0:lunpck), mpb(0:lunpck), mpc(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa), mwd(lmwa), mwe(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmpck/mpa, mpb, mpc
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmsums/mjsums
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
      COMMON /fmwa/mwd, mwe
! ..
! .. Data Statements ..
      DATA lchars/'+', '-', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', &
        '.', 'E', 'D', 'Q', 'M', 'e', 'd', 'q', 'm'/
      DATA ltypes/1, 1, 10*2, 3, 8*4/
      DATA lvals/1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9*0/
! ..
!             KW is the unit number for standard output from the
!                FM package.  This includes trace output and error
!                messages.

      kw = 6

!             MAXINT should be set to a very large integer, possibly
!                    the largest representable integer for the current
!                    machine.  For most 32-bit machines, MAXINT is set
!                    to  2**53 - 1 = 9.007D+15  when double precision
!                    arithmetic is used for M-variables.  Using integer
!                    M-variables usually gives MAXINT = 2**31 - 1 =
!                    2 147 483 647.

!                    Setting MAXINT to a smaller number is ok, but this
!                    unnecessarily restricts the permissible range of
!                    MBASE and MXEXP.

!                    The following code should set MAXINT to the largest
!                    representable number of the form 2**J - 1.

!             The FMMSET call keeps some compilers from doing the 110
!             loop at the highest precision available and then rounding
!             to the declared precision.

      maxint = 3
10    CALL fmmset(maxint,ml,mld2,mlm1)
      IF (mld2==maxint .AND. mlm1/=ml) THEN
        maxint = ml
        GO TO 10
      END IF

!             INTMAX is a large value close to the overflow threshold
!                    for integer variables.  It is usually 2**31 - 1
!                    for machines with 32-bit integer arithmetic.

!          WARNING:  This loop causes integer overflow to occur, so it
!                    is a likely place for the program to fail when
!                    run on a different machine.  The loop below has
!                    been used successfully with Fortran compilers
!                    for many different machines, but even different
!                    versions of the same compiler may give different
!                    results.  Check the values of MAXINT and INTMAX
!                    if there are problems installing FM.

      intmax = 3
20    l = 2*intmax + 1
      IF (int(l/2)==intmax) THEN
        intmax = l
        GO TO 20
      END IF

!             DPMAX should be set to a value near the machine's double
!                   precision overflow threshold, so that DPMAX and
!                   1.0D0/DPMAX are both representable in double
!                   precision.

      dpmax = 1.0D+74

!             SPMAX should be set to a value near the machine's single
!                   precision overflow threshold, so that 1.01*SPMAX
!                   and 1.0/SPMAX are both representable in single
!                   precision.

      spmax = 1.0E+37

!             NDG2MX is the maximum value for NDIG that can be used
!                    internally.  FM routines may raise NDIG above
!                    NDIGMX temporarily, to compute correctly
!                    rounded results.
!                    In the definition of LUNPCK, the '6/5' condition
!                    allows for converting from a large base to the
!                    (smaller) largest power of ten base for output
!                    conversion.
!                    The '+ 20' condition allows for the need to carry
!                    many guard digits when using a small base like 2.

      ndg2mx = lunpck - 1

!             MXBASE is the maximum value for MBASE.

      temp = maxint
      mxbase = int(min(dble(intmax),sqrt(temp)))

!             MBASE is the currently used base for arithmetic.

      k = int(log10(dble(mxbase)/4))
      mbase = 10**k

!             NDIG is the number of digits currently being carried.

      npsave = nprec
      ndig = 2 + (nprec+2)/k
      IF (ndig<2 .OR. ndig>ndigmx) THEN
        ndig = max(2,min(ndigmx,ndig))
        WRITE (kw,90000) nprec, ndig
        npsave = 0
      END IF

!             KFLAG is the flag for error conditions.

      kflag = 0

!             NTRACE is the trace switch.  Default is no printing.

      ntrace = 0

!             LVLTRC is the trace level.  Default is to trace only
!                    routines called directly by the user.

      lvltrc = 1

!             NCALL is the call stack pointer.

      ncall = 0

!             NAMEST is the call stack.

      DO 30 j = 0, 50
        namest(j) = 'MAIN  '
30    CONTINUE

!             Some constants that are often needed are stored with the
!             maximum precision to which they have been computed in the
!             currently used base.  This speeds up the trig, log, power,
!             and exponential functions.

!             NDIGPI is the number of digits available in the currently
!                    stored value of pi (MPISAV).

      ndigpi = 0

!             MBSPI is the value of MBASE for the currently stored
!                   value of pi.

      mbspi = 0

!             NDIGE is the number of digits available in the currently
!                   stored value of e (MESAV).

      ndige = 0

!             MBSE is the value of MBASE for the currently stored
!                  value of e.

      mbse = 0

!             NDIGLB is the number of digits available in the currently
!                    stored value of LN(MBASE) (MLBSAV).

      ndiglb = 0

!             MBSLB is the value of MBASE for the currently stored
!                   value of LN(MBASE).

      mbslb = 0

!             NDIGLI is the number of digits available in the currently
!                    stored values of the four logarithms used by FMLNI
!                    MLN1 - MLN4.

      ndigli = 0

!             MBSLI is the value of MBASE for the currently stored
!                   values of MLN1 - MLN4.

      mbsli = 0

!             MXEXP  is the current maximum exponent.
!             MXEXP2 is the internal maximum exponent. This is used to
!                    define the overflow and underflow thresholds.

!             These values are chosen so that FM routines can raise the
!             overflow/underflow limit temporarily while computing
!             intermediate results, and so that EXP(INTMAX) is greater
!             than MXBASE**(MXEXP2+1).

!             The overflow threshold is MBASE**(MXEXP+1), and the
!             underflow threshold is MBASE**(-MXEXP-1).
!             This means the valid exponents in the first word of an FM
!             number can range from -MXEXP to MXEXP+1 (inclusive).

      mxexp = int((dble(intmax))/(2.0D0*log(dble(mxbase)))-1.0D0)
      mxexp2 = int(2*mxexp+mxexp/100)

!             KACCSW is a switch used to enable cancellation error
!                    monitoring.  Routines where cancellation is
!                    not a problem run faster by skipping the
!                    cancellation monitor calculations.
!                    KACCSW = 0 means no error monitoring,
!                           = 1 means error monitoring is done.

      kaccsw = 0

!             MEXPUN is the exponent used as a special symbol for
!                    underflowed results.

      mexpun = -mxexp2 - 5*ndigmx

!             MEXPOV is the exponent used as a special symbol for
!                    overflowed results.

      mexpov = -mexpun

!             MUNKNO is the exponent used as a special symbol for
!                    unknown FM results (1/0, SQRT(-3.0), ...).

      munkno = mexpov + 5*ndigmx

!             RUNKNO is returned from FM to real or double conversion
!                    routines when no valid result can be expressed in
!                    real or double precision.  On systems that provide
!                    a value for undefined results (e.g., Not A Number)
!                    setting RUNKNO to that value is reasonable.  On
!                    other systems set it to a value that is likely to
!                    make any subsequent results obviously wrong that
!                    use it.  In either case a KFLAG = -4 condition is
!                    also returned.

      runkno = -1.01*spmax

!             IUNKNO is returned from FM to integer conversion routines
!                    when no valid result can be expressed as a one word
!                    integer.  KFLAG = -4 is also set.

      iunkno = -int(mxexp2)

!             JFORM1 indicates the format used by FMOUT.

      jform1 = 1

!             JFORM2 indicates the number of digits used in FMOUT.

      jform2 = npsave

!             KRAD = 1 indicates that trig functions use radians,
!                  = 0 means use degrees.

      krad = 1

!             KWARN = 0 indicates that no warning message is printed
!                       and execution continues when UNKNOWN or another
!                       exception is produced.
!                   = 1 means print a warning message and continue.
!                   = 2 means print a warning message and stop.

      kwarn = 1

!             KROUND = 1   causes all results to be rounded to the
!                          nearest FM number, or to the value with
!                          an even last digit if the result is halfway
!                          between two FM numbers.
!                    = 0   causes all results to be chopped.

      kround = 1

!             KSWIDE defines the maximum screen width to be used for
!                    all unit KW output.

      kswide = 80

!             KESWCH = 1   causes input to FMINP with no digits before
!                          the exponent letter to be treated as if there
!                          were a leading '1'.  This is sometimes better
!                          for interactive input:  'E7' converts to
!                          10.0**7.
!                    = 0   causes a leading zero to be assumed.  This
!                          gives compatibility with Fortran:  'E7'
!                          converts to 0.0.

      keswch = 1

!             CMCHAR defines the exponent letter to be used for
!                    FM variable output from FMOUT, as in 1.2345M+678.
!                    Change it to 'E' for output to be read by a
!                    non-FM program.

      cmchar = 'M'

!             KSUB is an internal flag set during subtraction so that
!                  the addition routine will negate its second argument.

      ksub = 0

!             KDEBUG = 0   Error checking is not done for valid input
!                          arguments and parameters like NDIG and MBASE
!                          upon entry to each routine.
!                    = 1   Error checking is done.

      kdebug = 0

!             Initialize two hash tables that are used for character
!             look-up during input conversion.

      DO 40 j = lhash1, lhash2
        khasht(j) = 5
        khashv(j) = 0
40    CONTINUE
      DO 50 j = 1, 21
        kpt = ichar(lchars(j))
        IF (kpt<lhash1 .OR. kpt>lhash2) THEN
          WRITE (kw,90010) lchars(j), kpt, lhash1, lhash2
        ELSE
          khasht(kpt) = ltypes(j)
          khashv(kpt) = lvals(j)
        END IF
50    CONTINUE

!             DPEPS is the approximate machine precision.

      one = 1.0D0
      two = 128.0D0
      dpeps = one

60    dpeps = dpeps/two
      CALL fmdbl(one,dpeps,yt)
      IF (yt>one) GO TO 60
      dpeps = dpeps*two
      two = 2.0D0
70    dpeps = dpeps/two
      CALL fmdbl(one,dpeps,yt)
      IF (yt>one) GO TO 70
      dpeps = dpeps*two

!             FMCONS sets several real and double precision constants.

      CALL fmcons

      RETURN
90000 FORMAT (//' Precision out of range when calling FMSET.','  NPREC =', &
        I20/' The nearest valid NDIG will be used',' instead:   NDIG =',I6//)
90010 FORMAT (/' Error in input conversion.'/ &
        ' ICHAR function was out of range for the current', &
        ' dimensions.'/' ICHAR(''',A,''') gave the value ',I12, &
        ', which is outside the currently'/' dimensioned',' bounds of (',I5, &
        ':',I5,') for variables KHASHT ','and KHASHV.'/ &
        ' Re-define the two parameters ', &
        'LHASH1 and LHASH2 so the dimensions will'/' contain', &
        ' all possible output values from ICHAR.'//)
    END SUBROUTINE fmset
    SUBROUTINE fmabs(ma,mb)

!  MB = ABS(MA)

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, log, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: md2b
      INTEGER :: kwrnsv
! ..
! .. External Subroutines ..
      EXTERNAL fmeq, fmntr
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      namest(ncall) = 'FMABS '
      IF (ntrace/=0) CALL fmntr(2,ma,ma,1)

      kflag = 0
      kwrnsv = kwarn
      kwarn = 0
      CALL fmeq(ma,mb)
      mb(2) = abs(mb(2))
      kwarn = kwrnsv

      IF (kaccsw==1) THEN
        md2b = nint((ndig-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
        mb(0) = min(mb(0),md2b)
      END IF
      IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmabs
    SUBROUTINE fmacos(ma,mb)

!  MB = ACOS(MA)

      IMPLICIT NONE

!             Scratch array usage during FMACOS:   M01 - M06

! .. Intrinsic Functions ..
      INTRINSIC abs, log, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, macca, macmax, mxsave
      INTEGER :: k, kasave, kovun, kreslt, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmatan, fmcons, fmdiv, fmentr, fmeq2, fmexit, fmi2m, &
        fmmpy, fmntr, fmpi, fmrslt, fmsqrt, fmsub, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (abs(ma(1))>mexpab .OR. ma(1)>0 .OR. ma(2)==0) THEN
        CALL fmentr('FMACOS',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        namest(ncall) = 'FMACOS'
        IF (ntrace/=0) CALL fmntr(2,ma,ma,1)
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
        ndsave = ndig
        IF (ncall==1) THEN
          k = max(ngrd52-1,2)
          ndig = max(ndig+k,2)
          IF (ndig>ndg2mx) THEN
            kflag = -9
            CALL fmwarn
            ndig = ndsave
            kreslt = 12
            CALL fmrslt(ma,ma,mb,kreslt)
            IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
            ncall = ncall - 1
            RETURN
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 0
        mxsave = mxexp
        mxexp = mxexp2
      END IF

      ma2 = ma(2)
      macca = ma(0)
      CALL fmeq2(ma,mb,ndsave,ndig,0)
      mb(0) = nint(ndig*alogm2)

!             Use ACOS(X) = ATAN(SQRT(1-X*X)/X)

      mb(2) = abs(mb(2))
      CALL fmi2m(1,m05)
      CALL fmsub(m05,mb,m03)
      CALL fmadd(m05,mb,m04)
      CALL fmmpy(m03,m04,m04)
      CALL fmsqrt(m04,m04)
      CALL fmdiv(m04,mb,mb)

      CALL fmatan(mb,mb)

      IF (ma2<0) THEN
        IF (krad==1) THEN
          CALL fmpi(m05)
        ELSE
          CALL fmi2m(180,m05)
        END IF
        CALL fmsub(m05,mb,mb)
      END IF

!             Round the result and return.

      macmax = nint((ndsave-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
      mb(0) = min(mb(0),macca,macmax)
      CALL fmexit(mb,mb,ndsave,mxsave,kasave,kovun)
      RETURN
    END SUBROUTINE fmacos
    SUBROUTINE fmadd(ma,mb,mc)

!  MC = MA + MB

!  This routine performs the trace printing for addition.
!  FMADD2 is used to do the arithmetic.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. External Subroutines ..
      EXTERNAL fmadd2, fmntr
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (ntrace/=0) THEN
        namest(ncall) = 'FMADD '
        CALL fmntr(2,ma,mb,2)

        CALL fmadd2(ma,mb,mc)

        CALL fmntr(1,mc,mc,1)
      ELSE
        CALL fmadd2(ma,mb,mc)
      END IF
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmadd
    SUBROUTINE fmadd2(ma,mb,mc)

!  Internal addition routine.  MC = MA + MB

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, log, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      REAL :: b2rda, b2rdb
      REAL (KIND(0.0D0)) :: ma0, ma1, ma2, mb0, mb1, mb2, mb2rd
      INTEGER :: j, jcomp, jsign, kreslt, n1, nguard, nmwa
! ..
! .. External Subroutines ..
      EXTERNAL fmaddn, fmaddp, fmargs, fmcons, fmeq, fmmove, fmrslt, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (abs(ma(1))>mexpab .OR. abs(mb(1))>mexpab .OR. kdebug==1) THEN
        IF (ksub==1) THEN
          CALL fmargs('FMSUB ',2,ma,mb,kreslt)
        ELSE
          CALL fmargs('FMADD ',2,ma,mb,kreslt)
        END IF
        IF (kreslt/=0) THEN
          ncall = ncall + 1
          IF (ksub==1) THEN
            namest(ncall) = 'FMSUB '
          ELSE
            namest(ncall) = 'FMADD '
          END IF
          CALL fmrslt(ma,mb,mc,kreslt)
          ncall = ncall - 1
          RETURN
        END IF
      ELSE
        IF (ma(2)==0) THEN
          ma0 = min(ma(0),mb(0))
          CALL fmeq(mb,mc)
          mc(0) = ma0
          kflag = 1
          IF (ksub==1) THEN
            IF (mc(1)/=munkno) mc(2) = -mc(2)
            kflag = 0
          END IF
          RETURN
        END IF
        IF (mb(2)==0) THEN
          ma0 = min(ma(0),mb(0))
          CALL fmeq(ma,mc)
          mc(0) = ma0
          kflag = 1
          RETURN
        END IF
      END IF

      ma0 = ma(0)
      IF (kaccsw==1) THEN
        mb0 = mb(0)
        ma1 = ma(1)
        mb1 = mb(1)
      END IF
      kflag = 0
      n1 = ndig + 1

!             NGUARD is the number of guard digits used.

      IF (ncall>1) THEN
        nguard = ngrd21
        IF (nguard>ndig) nguard = ndig
      ELSE
        nguard = ngrd52
        IF (nguard>ndig) nguard = ndig
      END IF
      nmwa = n1 + nguard

!             Save the signs of MA and MB and then work with
!             positive numbers.
!             JSIGN is the sign of the result of MA + MB.

      jsign = 1
      ma2 = ma(2)
      mb2 = mb(2)
      IF (ksub==1) mb2 = -mb2
      ma(2) = abs(ma(2))
      mb(2) = abs(mb(2))

!             See which one is larger in absolute value.

      IF (ma(1)>mb(1)) THEN
        jcomp = 1
        GO TO 20
      END IF
      IF (mb(1)>ma(1)) THEN
        jcomp = 3
        GO TO 20
      END IF

      DO 10 j = 2, n1
        IF (ma(j)>mb(j)) THEN
          jcomp = 1
          GO TO 20
        END IF
        IF (mb(j)>ma(j)) THEN
          jcomp = 3
          GO TO 20
        END IF
10    CONTINUE

      jcomp = 2

20    IF (jcomp<3) THEN
        IF (ma2<0) jsign = -1
        IF (ma2*mb2>0) THEN
          CALL fmaddp(ma,mb,nguard,nmwa)
        ELSE
          CALL fmaddn(ma,mb,nguard,nmwa)
        END IF
      ELSE
        IF (mb2<0) jsign = -1
        IF (ma2*mb2>0) THEN
          CALL fmaddp(mb,ma,nguard,nmwa)
        ELSE
          CALL fmaddn(mb,ma,nguard,nmwa)
        END IF
      END IF
      IF (ksub==1) mb2 = -mb2
      mb(2) = mb2
      ma(2) = ma2

!             Transfer to MC and fix the sign of the result.

      CALL fmmove(mwa,mc)
      IF (jsign<0) mc(2) = -mc(2)

      IF (kflag<0) THEN
        IF (ksub==1) THEN
          namest(ncall) = 'FMSUB '
        ELSE
          namest(ncall) = 'FMADD '
        END IF
        CALL fmwarn
      END IF

      IF (kaccsw==1) THEN
        b2rda = log(real(abs(mc(2))+1)/real(abs(ma2)+1))/0.69315 + &
          real(mc(1)-ma1)*alogm2 + real(ma0)
        b2rdb = log(real(abs(mc(2))+1)/real(abs(mb2)+1))/0.69315 + &
          real(mc(1)-mb1)*alogm2 + real(mb0)
        mb2rd = nint(max(0.0,min(b2rda,b2rdb,(ndig-1)*alogm2+log(real(abs(mc(2 &
          ))+1))/0.69315)))
        IF (mc(2)==0) THEN
          mc(0) = 0
        ELSE
          mc(0) = min(max(ma0,mb0),mb2rd)
        END IF
      ELSE
        mc(0) = ma0
      END IF

      RETURN
    END SUBROUTINE fmadd2
    SUBROUTINE fmaddi(ma,ival)

!  MA = MA + IVAL

!  Increment MA by one word integer IVAL.

!  This routine is faster than FMADD when IVAL is small enough so
!  that it can be added to a single word of MA without often causing
!  a carry.  Otherwise FMI2M and FMADD are used.

      IMPLICIT NONE

!             Scratch array usage during FMADDI:   M01

! .. Intrinsic Functions ..
      INTRINSIC abs, int, log, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maexp, md2b, mksum
      INTEGER :: kptma
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmi2m, fmntr, fmntri
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (ntrace/=0) THEN
        namest(ncall) = 'FMADDI'
        CALL fmntr(2,ma,ma,1)
        CALL fmntri(2,ival,0)
      END IF
      kflag = 0

      maexp = ma(1)
      IF (maexp<=0 .OR. maexp>ndig) GO TO 10
      kptma = int(maexp) + 1
      IF (kptma>2 .AND. ma(2)<0) THEN
        mksum = ma(kptma) - ival
      ELSE
        mksum = ma(kptma) + ival
      END IF

      IF (mksum>=mbase .OR. mksum<=(-mbase)) GO TO 10
      IF (ma(2)<0) THEN
        IF (kptma>2) THEN
          IF (mksum>=0) THEN
            ma(kptma) = mksum
            GO TO 20
          ELSE
            GO TO 10
          END IF
        ELSE
          IF (mksum<0) THEN
            ma(kptma) = mksum
            GO TO 20
          ELSE
            GO TO 10
          END IF
        END IF
      ELSE
        IF (kptma>2) THEN
          IF (mksum>=0) THEN
            ma(kptma) = mksum
            GO TO 20
          ELSE
            GO TO 10
          END IF
        ELSE
          IF (mksum>0) THEN
            ma(kptma) = mksum
            GO TO 20
          ELSE
            GO TO 10
          END IF
        END IF
      END IF

10    CALL fmi2m(ival,m01)
      CALL fmadd(ma,m01,ma)

20    IF (kaccsw==1) THEN
        md2b = nint((ndig-1)*alogm2+log(real(abs(ma(2))+1))/0.69315)
        ma(0) = min(ma(0),md2b)
      END IF
      IF (ntrace/=0) THEN
        CALL fmntr(1,ma,ma,1)
      END IF
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmaddi
    SUBROUTINE fmaddn(ma,mb,nguard,nmwa)

!  Internal addition routine.  MWA = MA - MB
!  The arguments are such that MA.GE.MB.GE.0.

!  NGUARD is the number of guard digits being carried.
!  NMWA is the number of words in MWA that will be used.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, int, min
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: nguard, nmwa
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: mk, mr
      INTEGER :: j, k, kl, kp1, kp2, kpt, ksh, n1, n2, nk, nk1
! ..
! .. External Subroutines ..
      EXTERNAL fmrnd
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      n1 = ndig + 1

!             Check for an insignificant operand.

      mk = ma(1) - mb(1)
      IF (mk>=ndig+2) THEN
        DO 10 j = 1, n1
          mwa(j) = ma(j)
10      CONTINUE
        mwa(n1+1) = 0
        kflag = 1
        RETURN
      END IF
      k = int(mk)
      IF (nguard<=1) nmwa = n1 + 2

!             Subtract MB from MA.

      kp1 = min(n1,k+1)
      mwa(k+1) = 0
      DO 20 j = 1, kp1
        mwa(j) = ma(j)
20    CONTINUE
      kp2 = k + 2

!             (Inner Loop)

      DO 30 j = kp2, n1
        mwa(j) = ma(j) - mb(j-k)
30    CONTINUE

      n2 = ndig + 2
      IF (n2-k<=1) n2 = 2 + k
      nk = min(nmwa,n1+k)
      DO 40 j = n2, nk
        mwa(j) = -mb(j-k)
40    CONTINUE
      nk1 = nk + 1
      DO 50 j = nk1, nmwa
        mwa(j) = 0
50    CONTINUE

!             Normalize.  Fix the sign of any negative digit.

      IF (k>0) THEN
        DO 60 j = nmwa, kp2, -1
          IF (mwa(j)<0) THEN
            mwa(j) = mwa(j) + mbase
            mwa(j-1) = mwa(j-1) - 1
          END IF
60      CONTINUE

        kpt = kp2 - 1
70      IF (mwa(kpt)<0 .AND. kpt>=3) THEN
          mwa(kpt) = mwa(kpt) + mbase
          mwa(kpt-1) = mwa(kpt-1) - 1
          kpt = kpt - 1
          GO TO 70
        END IF
        GO TO 90
      END IF

      DO 80 j = n1, 3, -1
        IF (mwa(j)<0) THEN
          mwa(j) = mwa(j) + mbase
          mwa(j-1) = mwa(j-1) - 1
        END IF
80    CONTINUE

!             Shift left if there are any leading zeros in the mantissa.

90    DO 100 j = 2, nmwa
        IF (mwa(j)>0) THEN
          ksh = j - 2
          GO TO 110
        END IF
100   CONTINUE
      mwa(1) = 0
      RETURN

110   IF (ksh>0) THEN
        kl = nmwa - ksh
        DO 120 j = 2, kl
          mwa(j) = mwa(j+ksh)
120     CONTINUE
        DO 130 j = kl + 1, nmwa
          mwa(j) = 0
130     CONTINUE
        mwa(1) = mwa(1) - ksh
      END IF

!             Round the result.

      mr = 2*mwa(ndig+2) + 1
      IF (mr>=mbase) THEN
        IF (mr-1>mbase .AND. mwa(n1)<mbase-1) THEN
          IF (kround/=0 .OR. ncall>1) THEN
            mwa(n1) = mwa(n1) + 1
            mwa(n1+1) = 0
          END IF
        ELSE
          CALL fmrnd(mwa,ndig,nguard,0)
        END IF
      END IF

!             See if the result is equal to one of the input arguments.

      IF (abs(ma(1)-mb(1))<ndig) GO TO 150
      IF (abs(ma(1)-mb(1))>ndig+1) THEN
        kflag = 1
        GO TO 150
      END IF

      n2 = ndig + 4
      DO 140 j = 3, n1
        IF (mwa(n2-j)/=ma(n2-j)) GO TO 150
140   CONTINUE
      IF (mwa(1)/=ma(1)) GO TO 150
      IF (mwa(2)/=abs(ma(2))) GO TO 150
      kflag = 1

150   RETURN
    END SUBROUTINE fmaddn
    SUBROUTINE fmaddp(ma,mb,nguard,nmwa)

!  Internal addition routine.  MWA = MA + MB
!  The arguments are such that MA.GE.MB.GE.0.

!  NMWA is the number of words in MWA that will be used.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dint, int, min
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: nguard, nmwa
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: mk, mkt, mr
      INTEGER :: j, k, kp, kp2, kpt, kshift, n1, n2, nk
! ..
! .. External Subroutines ..
      EXTERNAL fmrnd
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      n1 = ndig + 1

!             Check for an insignificant operand.

      mk = ma(1) - mb(1)
      IF (mk>=ndig+1) THEN
        mwa(1) = ma(1) + 1
        mwa(2) = 0
        DO 10 j = 2, n1
          mwa(j+1) = ma(j)
10      CONTINUE
        mwa(n1+2) = 0
        kflag = 1
        RETURN
      END IF
      k = int(mk)

!             Add MA and MB.

      mwa(1) = ma(1) + 1
      mwa(2) = 0
      DO 20 j = 2, k + 1
        mwa(j+1) = ma(j)
20    CONTINUE
      kp2 = k + 2

!             (Inner Loop)

      DO 30 j = kp2, n1
        mwa(j+1) = ma(j) + mb(j-k)
30    CONTINUE
      n2 = ndig + 2
      nk = min(nmwa,n1+k)
      DO 40 j = n2, nk
        mwa(j+1) = mb(j-k)
40    CONTINUE
      DO 50 j = nk + 1, nmwa
        mwa(j+1) = 0
50    CONTINUE

!             Normalize.  Fix any digit not less than MBASE.

      IF (k==ndig) GO TO 120

      IF (k>0) THEN
        DO 60 j = n1 + 1, kp2, -1
          IF (mwa(j)>=mbase) THEN
            mwa(j) = mwa(j) - mbase
            mwa(j-1) = mwa(j-1) + 1
          END IF
60      CONTINUE

        kpt = kp2 - 1
70      IF (mwa(kpt)>=mbase .AND. kpt>=3) THEN
          mwa(kpt) = mwa(kpt) - mbase
          mwa(kpt-1) = mwa(kpt-1) + 1
          kpt = kpt - 1
          GO TO 70
        END IF
        GO TO 90
      END IF

      DO 80 j = n1 + 1, 3, -1
        IF (mwa(j)>=mbase) THEN
          mwa(j) = mwa(j) - mbase
          mwa(j-1) = mwa(j-1) + 1
        END IF
80    CONTINUE

!             Shift right if the leading digit is not less than MBASE.

90    IF (mwa(2)>=mbase) THEN
100     kp = nmwa + 4
        DO 110 j = 4, nmwa
          mwa(kp-j) = mwa(kp-j-1)
110     CONTINUE
        mkt = dint(mwa(2)/mbase)
        mwa(3) = mwa(2) - mkt*mbase
        mwa(2) = mkt
        mwa(1) = mwa(1) + 1
        IF (mwa(2)>=mbase) GO TO 100
      END IF

!             Round the result.

120   kshift = 0
      IF (mwa(2)==0) kshift = 1
      mr = 2*mwa(ndig+2+kshift) + 1
      IF (mr>=mbase) THEN
        IF (mr-1>mbase .AND. mwa(n1+kshift)<mbase-1) THEN
          IF (kround/=0 .OR. ncall>1) THEN
            mwa(n1+kshift) = mwa(n1+kshift) + 1
            mwa(n1+1+kshift) = 0
          END IF
        ELSE
          CALL fmrnd(mwa,ndig,nguard,kshift)
        END IF
      END IF

!             See if the result is equal to one of the input arguments.

      IF (abs(ma(1)-mb(1))<ndig) GO TO 140
      IF (kshift==0) GO TO 140
      IF (abs(ma(1)-mb(1))>ndig+1) THEN
        kflag = 1
        GO TO 140
      END IF

      n2 = ndig + 4
      DO 130 j = 3, n1
        IF (mwa(n2-j+1)/=ma(n2-j)) GO TO 140
130   CONTINUE
      IF (mwa(1)/=ma(1)+1) GO TO 140
      IF (mwa(3)/=abs(ma(2))) GO TO 140
      kflag = 1

140   RETURN
    END SUBROUTINE fmaddp
    SUBROUTINE fmargs(kroutn,nargs,ma,mb,kreslt)

!  Check the input arguments to a routine for special cases.

!  KROUTN - Name of the subroutine that was called
!  NARGS  - The number of input arguments (1 or 2)
!  MA     - First input argument
!  MB     - Second input argument (if NARGS is 2)
!  KRESLT - Result code returned to the calling routine.

!  Result codes:

!   0 - Perform the normal operation
!   1 - The result is the first input argument
!   2 - The result is the second input argument
!   3 - The result is -OVERFLOW
!   4 - The result is +OVERFLOW
!   5 - The result is -UNDERFLOW
!   6 - The result is +UNDERFLOW
!   7 - The result is -1.0
!   8 - The result is +1.0
!   9 - The result is -pi/2
!  10 - The result is +pi/2
!  11 - The result is 0.0
!  12 - The result is UNKNOWN
!  13 - The result is +pi
!  14 - The result is -pi/4
!  15 - The result is +pi/4

      IMPLICIT NONE

!             These tables define the result codes to be returned for
!             given values of the input argument(s).

!             For example, in row 7 column 2 of this DATA statement
!             KADD(2,7) = 2 means that if the first argument in a call
!             to FMADD is in category 7 ( -UNDERFLOW ) and the second
!             argument is in category 2 ( near -OVERFLOW but
!             representable ) then the result code is 2 ( the value
!             of the sum is equal to the second input argument).
!             See routine FMCAT for descriptions of the categories.

! .. Intrinsic Functions ..
      INTRINSIC abs, int, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kreslt, nargs
      CHARACTER (6) :: kroutn
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: mbs
      INTEGER :: j, kwrnsv, ncatma, ncatmb, nds
! ..
! .. Local Arrays ..
      INTEGER :: kacos(15), kadd(15,15), kasin(15), katan(15), kcos(15), &
        kcosh(15), kdiv(15,15), kexp(15), klg10(15), kln(15), kmpy(15,15), &
        kpwr(15,15), ksin(15), ksinh(15), ksqrt(15), ktan(15), ktanh(15)
! ..
! .. External Subroutines ..
      EXTERNAL fmcat, fmcons, fmim, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
! .. Data Statements ..
      DATA kadd/3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 12, 12, 3, 0, 0, 0, 0, &
        0, 1, 1, 1, 0, 0, 0, 0, 0, 12, 3, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, &
        0, 4, 3, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 4, 3, 0, 0, 0, 0, 0, &
        1, 1, 1, 0, 0, 0, 0, 0, 4, 3, 0, 0, 0, 0, 0, 12, 1, 12, 0, 0, 0, 0, 0, &
        4, 3, 2, 2, 2, 2, 12, 12, 5, 12, 12, 2, 2, 2, 2, 4, 3, 2, 2, 2, 2, 2, &
        5, 2, 6, 2, 2, 2, 2, 2, 4, 3, 2, 2, 2, 2, 12, 12, 6, 12, 12, 2, 2, 2, &
        2, 4, 3, 0, 0, 0, 0, 0, 12, 1, 12, 0, 0, 0, 0, 0, 4, 3, 0, 0, 0, 0, 0, &
        1, 1, 1, 0, 0, 0, 0, 0, 4, 3, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, &
        4, 3, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 4, 12, 0, 0, 0, 0, 0, 1, &
        1, 1, 0, 0, 0, 0, 0, 4, 12, 12, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4/
      DATA kmpy/4, 4, 4, 4, 12, 12, 12, 11, 12, 12, 12, 3, 3, 3, 3, 4, 0, 0, &
        0, 0, 0, 12, 11, 12, 0, 0, 1, 0, 0, 3, 4, 0, 0, 0, 0, 0, 12, 11, 12, &
        0, 0, 1, 0, 0, 3, 4, 0, 0, 0, 0, 0, 6, 11, 5, 0, 0, 1, 0, 0, 3, 12, 0, &
        0, 0, 0, 0, 6, 11, 5, 0, 0, 1, 0, 0, 12, 12, 0, 0, 0, 0, 0, 6, 11, 5, &
        0, 0, 1, 0, 0, 12, 12, 12, 12, 6, 6, 6, 6, 11, 5, 5, 5, 5, 12, 12, 12, &
        11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, &
        12, 5, 5, 5, 5, 11, 6, 6, 6, 6, 12, 12, 12, 12, 0, 0, 0, 0, 0, 5, 11, &
        6, 0, 0, 1, 0, 0, 12, 12, 0, 0, 0, 0, 0, 5, 11, 6, 0, 0, 1, 0, 0, 12, &
        3, 2, 2, 2, 2, 2, 5, 11, 6, 2, 2, 2, 2, 2, 4, 3, 0, 0, 0, 0, 0, 12, &
        11, 12, 0, 0, 1, 0, 0, 4, 3, 0, 0, 0, 0, 0, 12, 11, 12, 0, 0, 1, 0, 0, &
        4, 3, 3, 3, 3, 12, 12, 12, 11, 12, 12, 12, 4, 4, 4, 4/
      DATA kdiv/12, 12, 12, 4, 4, 4, 4, 12, 3, 3, 3, 3, 12, 12, 12, 12, 0, 0, &
        0, 0, 0, 4, 12, 3, 0, 0, 1, 0, 0, 12, 12, 0, 0, 0, 0, 0, 4, 12, 3, 0, &
        0, 1, 0, 0, 12, 6, 0, 0, 0, 0, 0, 4, 12, 3, 0, 0, 1, 0, 0, 5, 6, 0, 0, &
        0, 0, 0, 12, 12, 12, 0, 0, 1, 0, 0, 5, 6, 0, 0, 0, 0, 0, 12, 12, 12, &
        0, 0, 1, 0, 0, 5, 6, 6, 6, 6, 12, 12, 12, 12, 12, 12, 12, 5, 5, 5, 5, &
        11, 11, 11, 11, 11, 11, 11, 12, 11, 11, 11, 11, 11, 11, 11, 5, 5, 5, &
        5, 12, 12, 12, 12, 12, 12, 12, 6, 6, 6, 6, 5, 0, 0, 0, 0, 0, 12, 12, &
        12, 0, 0, 1, 0, 0, 6, 5, 0, 0, 0, 0, 0, 12, 12, 12, 0, 0, 1, 0, 0, 6, &
        5, 0, 0, 0, 0, 0, 3, 12, 4, 0, 0, 1, 0, 0, 6, 12, 0, 0, 0, 0, 0, 3, &
        12, 4, 0, 0, 1, 0, 0, 12, 12, 0, 0, 0, 0, 0, 3, 12, 4, 0, 0, 1, 0, 0, &
        12, 12, 12, 12, 3, 3, 3, 3, 12, 4, 4, 4, 4, 12, 12, 12/
      DATA kpwr/12, 12, 0, 5, 12, 12, 12, 8, 12, 12, 12, 3, 0, 12, 12, 12, 12, &
        0, 0, 12, 12, 12, 8, 12, 12, 12, 1, 0, 12, 12, 12, 12, 0, 0, 12, 12, &
        12, 8, 12, 12, 12, 1, 0, 12, 12, 12, 12, 0, 0, 12, 12, 12, 8, 12, 12, &
        12, 1, 0, 12, 12, 12, 12, 0, 0, 12, 12, 12, 8, 12, 12, 12, 1, 0, 12, &
        12, 12, 12, 0, 0, 12, 12, 12, 8, 12, 12, 12, 1, 0, 12, 12, 12, 12, 0, &
        3, 12, 12, 12, 8, 12, 12, 12, 5, 0, 12, 12, 12, 12, 12, 12, 12, 12, &
        12, 12, 11, 11, 11, 11, 11, 11, 11, 4, 4, 4, 4, 12, 12, 12, 8, 12, 12, &
        12, 6, 6, 6, 6, 4, 4, 0, 0, 0, 8, 8, 8, 8, 0, 0, 1, 0, 6, 6, 4, 4, 0, &
        0, 0, 8, 8, 8, 8, 0, 0, 1, 0, 6, 6, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, &
        8, 8, 8, 8, 6, 6, 0, 0, 0, 8, 8, 8, 8, 8, 0, 1, 0, 4, 4, 6, 6, 0, 0, &
        0, 8, 8, 8, 8, 8, 0, 1, 0, 4, 4, 6, 6, 6, 6, 12, 12, 12, 8, 12, 12, &
        12, 4, 4, 4, 4/
      DATA ksqrt/12, 12, 12, 12, 12, 12, 12, 11, 12, 0, 0, 8, 0, 0, 12/
      DATA kexp/6, 6, 0, 0, 0, 8, 8, 8, 8, 8, 0, 0, 0, 4, 4/
      DATA kln/12, 12, 12, 12, 12, 12, 12, 12, 12, 0, 0, 11, 0, 0, 12/
      DATA ksin/12, 12, 0, 0, 0, 0, 5, 11, 6, 0, 0, 0, 0, 12, 12/
      DATA kcos/12, 12, 0, 0, 0, 8, 8, 8, 8, 8, 0, 0, 0, 12, 12/
      DATA ktan/12, 12, 0, 0, 0, 0, 5, 11, 6, 0, 0, 0, 0, 12, 12/
      DATA kasin/12, 12, 12, 9, 0, 0, 5, 11, 6, 0, 0, 10, 12, 12, 12/
      DATA kacos/12, 12, 12, 13, 0, 10, 10, 10, 10, 10, 0, 11, 12, 12, 12/
      DATA katan/9, 9, 0, 14, 0, 0, 5, 11, 6, 0, 0, 15, 0, 10, 10/
      DATA ksinh/3, 3, 0, 0, 0, 1, 5, 11, 6, 1, 0, 0, 0, 4, 4/
      DATA kcosh/4, 4, 0, 0, 0, 8, 8, 8, 8, 8, 0, 0, 0, 4, 4/
      DATA ktanh/7, 7, 0, 0, 0, 1, 5, 11, 6, 1, 0, 0, 0, 8, 8/
      DATA klg10/12, 12, 12, 12, 12, 12, 12, 12, 12, 0, 0, 11, 0, 0, 12/
! ..
      kreslt = 12
      kflag = -4
      IF (ma(1)==munkno) RETURN
      IF (nargs==2) THEN
        IF (mb(1)==munkno) RETURN
      END IF
      IF (mblogs/=mbase) CALL fmcons
      kflag = 0
      namest(ncall) = kroutn

!             Check the validity of parameters if this is a user call.

      IF (ncall>1 .AND. kdebug==0) GO TO 50

!             Check NDIG.

      IF (ndig<2 .OR. ndig>ndigmx) THEN
        kflag = -1
        CALL fmwarn
        nds = ndig
        IF (ndig<2) ndig = 2
        IF (ndig>ndigmx) ndig = ndigmx
        WRITE (kw,90000) nds, ndig
        RETURN
      END IF

!             Check MBASE.

      IF (mbase<2 .OR. mbase>mxbase) THEN
        kflag = -2
        CALL fmwarn
        mbs = mbase
        IF (mbase<2) mbase = 2
        IF (mbase>mxbase) mbase = mxbase
        WRITE (kw,90010) int(mbs), int(mbase)
        CALL fmcons
        RETURN
      END IF

!             Check exponent range.

      IF (ma(1)>mxexp+1 .OR. ma(1)<-mxexp) THEN
        IF (abs(ma(1))/=mexpov .OR. abs(ma(2))/=1) THEN
          CALL fmim(0,ma)
          kflag = -3
          CALL fmwarn
          ma(1) = munkno
          ma(2) = 1
          ma(0) = nint(ndig*alogm2)
          RETURN
        END IF
      END IF
      IF (nargs==2) THEN
        IF (mb(1)>mxexp+1 .OR. mb(1)<-mxexp) THEN
          IF (abs(mb(1))/=mexpov .OR. abs(mb(2))/=1) THEN
            CALL fmim(0,mb)
            kflag = -3
            CALL fmwarn
            mb(1) = munkno
            mb(2) = 1
            mb(0) = nint(ndig*alogm2)
            RETURN
          END IF
        END IF
      END IF

!             Check for properly normalized digits in the
!             input arguments.

      IF (abs(ma(1)-int(ma(1)))/=0) kflag = 1
      IF (ma(2)<=(-mbase) .OR. ma(2)>=mbase .OR. abs(ma(2)-int(ma(2)))/=0) &
        kflag = 2
      IF (kdebug==0) GO TO 20
      DO 10 j = 3, ndig + 1
        IF (ma(j)<0 .OR. ma(j)>=mbase .OR. abs(ma(j)-int(ma(j)))/=0) THEN
          kflag = j
          GO TO 20
        END IF
10    CONTINUE
20    IF (kflag/=0) THEN
        j = kflag
        mbs = ma(j)
        CALL fmim(0,ma)
        kflag = -4
        kwrnsv = kwarn
        IF (kwarn>=2) kwarn = 1
        CALL fmwarn
        kwarn = kwrnsv
        IF (kwarn>=1) THEN
          WRITE (kw,*) ' First invalid array element:  MA(', j, ') = ', mbs
        END IF
        ma(1) = munkno
        ma(2) = 1
        ma(0) = nint(ndig*alogm2)
        IF (kwarn>=2) THEN
          STOP
        END IF
        RETURN
      END IF
      IF (nargs==2) THEN
        IF (abs(mb(1)-int(mb(1)))/=0) kflag = 1
        IF (mb(2)<=(-mbase) .OR. mb(2)>=mbase .OR. abs(mb(2)-int(mb(2)))/=0) &
          kflag = 2
        IF (kdebug==0) GO TO 40
        DO 30 j = 3, ndig + 1
          IF (mb(j)<0 .OR. mb(j)>=mbase .OR. abs(mb(j)-int(mb(j)))/=0) THEN
            kflag = j
            GO TO 40
          END IF
30      CONTINUE
40      IF (kflag/=0) THEN
          j = kflag
          mbs = mb(j)
          CALL fmim(0,mb)
          kflag = -4
          kwrnsv = kwarn
          IF (kwarn>=2) kwarn = 1
          CALL fmwarn
          kwarn = kwrnsv
          IF (kwarn>=1) THEN
            WRITE (kw,*) ' First invalid array element:  MB(', j, ') = ', mbs
          END IF
          mb(1) = munkno
          mb(2) = 1
          mb(0) = nint(ndig*alogm2)
          IF (kwarn>=2) THEN
            STOP
          END IF
          RETURN
        END IF
      END IF

!             Check for special cases.

50    CALL fmcat(ma,ncatma)
      ncatmb = 0
      IF (nargs==2) CALL fmcat(mb,ncatmb)

      IF (kroutn=='FMADD ') THEN
        kreslt = kadd(ncatmb,ncatma)
        GO TO 60
      END IF

      IF (kroutn=='FMSUB ') THEN
        IF (ncatmb<16) ncatmb = 16 - ncatmb
        kreslt = kadd(ncatmb,ncatma)
        GO TO 60
      END IF

      IF (kroutn=='FMMPY ') THEN
        kreslt = kmpy(ncatmb,ncatma)
        GO TO 60
      END IF

      IF (kroutn=='FMDIV ') THEN
        kreslt = kdiv(ncatmb,ncatma)
        GO TO 60
      END IF

      IF (kroutn=='FMPWR ') THEN
        kreslt = kpwr(ncatmb,ncatma)
        GO TO 60
      END IF

      IF (kroutn=='FMSQRT') THEN
        kreslt = ksqrt(ncatma)
        GO TO 60
      END IF

      IF (kroutn=='FMEXP ') THEN
        kreslt = kexp(ncatma)
        GO TO 60
      END IF

      IF (kroutn=='FMLN  ') THEN
        kreslt = kln(ncatma)
        GO TO 60
      END IF

      IF (kroutn=='FMSIN ') THEN
        kreslt = ksin(ncatma)
        GO TO 60
      END IF

      IF (kroutn=='FMCOS ') THEN
        kreslt = kcos(ncatma)
        GO TO 60
      END IF

      IF (kroutn=='FMTAN ') THEN
        kreslt = ktan(ncatma)
        GO TO 60
      END IF

      IF (kroutn=='FMASIN') THEN
        kreslt = kasin(ncatma)
        IF ((ncatma==7 .OR. ncatma==9) .AND. krad==0) kreslt = 12
        GO TO 60
      END IF

      IF (kroutn=='FMACOS') THEN
        kreslt = kacos(ncatma)
        GO TO 60
      END IF

      IF (kroutn=='FMATAN') THEN
        kreslt = katan(ncatma)
        IF ((ncatma==7 .OR. ncatma==9) .AND. krad==0) kreslt = 12
        GO TO 60
      END IF

      IF (kroutn=='FMSINH') THEN
        kreslt = ksinh(ncatma)
        GO TO 60
      END IF

      IF (kroutn=='FMCOSH') THEN
        kreslt = kcosh(ncatma)
        GO TO 60
      END IF

      IF (kroutn=='FMTANH') THEN
        kreslt = ktanh(ncatma)
        GO TO 60
      END IF

      IF (kroutn=='FMLG10') THEN
        kreslt = klg10(ncatma)
        GO TO 60
      END IF

      kreslt = 0
      RETURN

60    IF (kreslt==12) THEN
        kflag = -4
        CALL fmwarn
      END IF
      IF (kreslt==3 .OR. kreslt==4) THEN
        IF (ncatma==1 .OR. ncatma==7 .OR. ncatma==9 .OR. ncatma==15 .OR. &
            ncatmb==1 .OR. ncatmb==7 .OR. ncatmb==9 .OR. ncatmb==15) THEN
          kflag = -5
        ELSE
          kflag = -5
          CALL fmwarn
        END IF
      END IF
      IF (kreslt==5 .OR. kreslt==6) THEN
        IF (ncatma==1 .OR. ncatma==7 .OR. ncatma==9 .OR. ncatma==15 .OR. &
            ncatmb==1 .OR. ncatmb==7 .OR. ncatmb==9 .OR. ncatmb==15) THEN
          kflag = -6
        ELSE
          kflag = -6
          CALL fmwarn
        END IF
      END IF
      RETURN
90000 FORMAT (' NDIG was',I10,'.  It has been changed to',I10,'.')
90010 FORMAT (' MBASE was',I10,'.  It has been changed to',I10,'.')
    END SUBROUTINE fmargs
    SUBROUTINE fmasin(ma,mb)

!  MB = ARCSIN(MA)

      IMPLICIT NONE

!             Scratch array usage during FMASIN:   M01 - M06

! .. Intrinsic Functions ..
      INTRINSIC abs, log, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: macca, macmax, mxsave
      INTEGER :: k, kasave, kovun, kreslt, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmatan, fmcons, fmdiv, fmentr, fmeq2, fmexit, fmi2m, &
        fmmpy, fmntr, fmrslt, fmsqrt, fmsub, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (abs(ma(1))>mexpab .OR. ma(1)>0 .OR. ma(2)==0) THEN
        CALL fmentr('FMASIN',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        namest(ncall) = 'FMASIN'
        IF (ntrace/=0) CALL fmntr(2,ma,ma,1)
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
        ndsave = ndig
        IF (ncall==1) THEN
          k = max(ngrd52-1,2)
          ndig = max(ndig+k,2)
          IF (ndig>ndg2mx) THEN
            kflag = -9
            CALL fmwarn
            ndig = ndsave
            kreslt = 12
            CALL fmrslt(ma,ma,mb,kreslt)
            IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
            ncall = ncall - 1
            RETURN
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 0
        mxsave = mxexp
        mxexp = mxexp2
      END IF

      macca = ma(0)
      CALL fmeq2(ma,mb,ndsave,ndig,0)
      mb(0) = nint(ndig*alogm2)

!             Use ASIN(X) = ATAN(X/SQRT(1-X*X))

      CALL fmi2m(1,m05)
      CALL fmsub(m05,mb,m03)
      CALL fmadd(m05,mb,m04)
      CALL fmmpy(m03,m04,m04)
      CALL fmsqrt(m04,m04)
      CALL fmdiv(mb,m04,mb)

      CALL fmatan(mb,mb)

!             Round the result and return.

      macmax = nint((ndsave-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
      mb(0) = min(mb(0),macca,macmax)
      CALL fmexit(mb,mb,ndsave,mxsave,kasave,kovun)
      RETURN
    END SUBROUTINE fmasin
    SUBROUTINE fmatan(ma,mb)

!  MB = ARCTAN(MA)

      IMPLICIT NONE

!             Scratch array usage during FMATAN:   M01 - M06

! .. Intrinsic Functions ..
      INTRINSIC abs, atan, int, log, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma1, ma2, macca, macmax, mxsave
      REAL (KIND(0.0D0)) :: x, xm
      INTEGER :: j, k, kasave, kovun, kreslt, krsave, kst, kwrnsv, ndsav1, &
        ndsave, ndsv
! ..
! .. Local Arrays ..
      INTEGER :: nstack(19)
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmdig, fmdiv, fmdivi, fmdpm, fmentr, fmeq, &
        fmeq2, fmexit, fmi2m, fmm2dp, fmmpy, fmmpyi, fmntr, fmpi, fmrslt, &
        fmsin, fmsqr, fmsqrt, fmsub, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (abs(ma(1))>mexpab .OR. ma(2)==0) THEN
        CALL fmentr('FMATAN',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        namest(ncall) = 'FMATAN'
        IF (ntrace/=0) CALL fmntr(2,ma,ma,1)
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
        ndsave = ndig
        IF (ncall==1) THEN
          k = max(ngrd52-1,2)
          ndig = max(ndig+k,2)
          IF (ndig>ndg2mx) THEN
            kflag = -9
            CALL fmwarn
            ndig = ndsave
            kreslt = 12
            CALL fmrslt(ma,ma,mb,kreslt)
            IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
            ncall = ncall - 1
            RETURN
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 0
        mxsave = mxexp
        mxexp = mxexp2
      END IF

      macca = ma(0)
      CALL fmeq2(ma,m05,ndsave,ndig,0)
      m05(0) = nint(ndig*alogm2)

!             If MA.GE.1 work with 1/MA.

      ma1 = ma(1)
      ma2 = ma(2)
      m05(2) = abs(m05(2))
      IF (ma1>=1) THEN
        CALL fmi2m(1,mb)
        CALL fmdiv(mb,m05,m05)
      END IF

      krsave = krad
      krad = 1
      kwrnsv = kwarn

      x = m05(1)
      xm = mxbase

!             In case pi has not been computed at the current precision
!             and will be needed here, get it to full precision first
!             to avoid repeated calls at increasing precision during
!             Newton iteration.

      IF (ma1>=1 .OR. krsave==0) THEN
        IF (mbspi/=mbase .OR. ndigpi<ndig) THEN
          ndsv = ndig
          ndig = min(ndig+2,ndg2mx)
          ncall = ncall + 1
          namest(ncall) = 'NOEQ  '
          CALL fmpi(mpisav)
          ncall = ncall - 1
          ndig = ndsv
        END IF
      END IF

!             If the argument is small, use the Taylor series,
!             otherwise use Newton iteration.

      IF (x*dlogmb<-5.0D0*log(xm)) THEN
        kwarn = 0
        CALL fmeq(m05,mb)
        IF (mb(1)<=-ndig) GO TO 30
        CALL fmsqr(m05,m06)
        j = 3
        ndsav1 = ndig

10      CALL fmmpy(m05,m06,m05)
        IF (m05(1)/=munkno) m05(2) = -m05(2)
        CALL fmdivi(m05,j,m03)
        ndig = ndsav1
        CALL fmadd(mb,m03,mb)
        IF (kflag/=0) THEN
          kflag = 0
          GO TO 30
        END IF
        ndig = ndsav1 - int((mb(1)-m03(1)))
        IF (ndig<2) ndig = 2
        j = j + 2
        GO TO 10
      ELSE

        CALL fmm2dp(m05,x)
        x = atan(x)
        CALL fmdpm(x,mb)
        CALL fmdig(nstack,kst)

!             Newton iteration.

        DO 20 j = 1, kst
          ndig = nstack(j)
          CALL fmsin(mb,m06)
          CALL fmsqr(m06,m03)
          CALL fmi2m(1,m04)
          CALL fmsub(m04,m03,m03)
          CALL fmsqrt(m03,m04)
          CALL fmdiv(m06,m04,m04)
          CALL fmsub(m04,m05,m04)
          CALL fmmpy(m03,m04,m04)
          CALL fmsub(mb,m04,mb)
20      CONTINUE
        mb(0) = nint(ndig*alogm2)
      END IF

!             If MA.GE.1 use pi/2 - ATAN(1/MA)

30    IF (ma1>=1) THEN
        CALL fmdivi(mpisav,2,m06)
        CALL fmsub(m06,mb,mb)
      END IF

!             Convert to degrees if necessary, round and return.

      krad = krsave
      IF (krad==0) THEN
        CALL fmmpyi(mb,180,mb)
        CALL fmdiv(mb,mpisav,mb)
      END IF
      IF (mb(1)/=munkno .AND. ma2<0) mb(2) = -mb(2)

      IF (kflag==1) kflag = 0
      kwarn = kwrnsv
      macmax = nint((ndsave-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
      mb(0) = min(mb(0),macca,macmax)
      CALL fmexit(mb,mb,ndsave,mxsave,kasave,kovun)
      RETURN
    END SUBROUTINE fmatan
    SUBROUTINE fmatn2(ma,mb,mc)

!  MC = ATAN2(MA,MB)

!  MC is returned as the angle between -pi and pi (or -180 and 180 if
!  degree mode is selected) for which TAN(MC) = MA/MB.  MC is an angle
!  for the point (MB,MA) in polar coordinates.

      IMPLICIT NONE

!             Scratch array usage during FMATN2:   M01 - M06

! .. Intrinsic Functions ..
      INTRINSIC abs, int, log, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: macca, maccb, macmax, mxexp1, mxsave
      INTEGER :: jquad, k, kasave, kovun, kreslt, kwrnsv, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmatan, fmcons, fmdiv, fmdivi, fmentr, fmeq2, fmexit, fmi2m, &
        fmim, fmntr, fmpi, fmrslt, fmsub, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (abs(ma(1))>mexpab .OR. abs(mb(1))>mexpab) THEN
        CALL fmentr('FMATN2',ma,mb,2,mc,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        namest(ncall) = 'FMATN2'
        IF (ntrace/=0) CALL fmntr(2,ma,mb,2)
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
        ndsave = ndig
        IF (ncall==1) THEN
          k = max(ngrd52-1,2)
          ndig = max(ndig+k,2)
          IF (ndig>ndg2mx) THEN
            kflag = -9
            CALL fmwarn
            ndig = ndsave
            kreslt = 12
            CALL fmrslt(ma,mb,mc,kreslt)
            IF (ntrace/=0) CALL fmntr(1,mc,mc,1)
            ncall = ncall - 1
            RETURN
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 0
        mxsave = mxexp
        mxexp = mxexp2
      END IF

      kwrnsv = kwarn
      kwarn = 0

      macca = ma(0)
      maccb = mb(0)
      CALL fmeq2(ma,m01,ndsave,ndig,0)
      m01(0) = nint(ndig*alogm2)
      CALL fmeq2(mb,m02,ndsave,ndig,0)
      m02(0) = nint(ndig*alogm2)

!             Check for special cases.

      IF (ma(1)==munkno .OR. mb(1)==munkno .OR. (ma(2)==0 .AND. mb(2)==0)) &
          THEN
        CALL fmim(0,mc)
        mc(1) = munkno
        mc(2) = 1
        mc(0) = nint(ndig*alogm2)
        kflag = -4
        GO TO 10
      END IF

      IF (mb(2)==0 .AND. ma(2)>0) THEN
        IF (krad==0) THEN
          CALL fmi2m(90,mc)
        ELSE
          CALL fmpi(mc)
          CALL fmdivi(mc,2,mc)
        END IF
        GO TO 10
      END IF

      IF (mb(2)==0 .AND. ma(2)<0) THEN
        IF (krad==0) THEN
          CALL fmi2m(-90,mc)
        ELSE
          CALL fmpi(mc)
          CALL fmdivi(mc,-2,mc)
        END IF
        GO TO 10
      END IF

      mxexp1 = int(mxexp2/2.01D0)
      IF (ma(1)==mexpov .AND. mb(1)<mxexp1-ndig-2) THEN
        IF (krad==0) THEN
          CALL fmi2m(90,mc)
        ELSE
          CALL fmpi(mc)
          CALL fmdivi(mc,2,mc)
        END IF
        IF (m01(2)<0) mc(2) = -mc(2)
        GO TO 10
      END IF

      IF (ma(1)==mexpun .AND. (-mb(1))<mxexp1-ndig-2 .AND. mb(2)<0) THEN
        IF (krad==0) THEN
          CALL fmi2m(180,mc)
        ELSE
          CALL fmpi(mc)
        END IF
        IF (m01(2)<0) mc(2) = -mc(2)
        GO TO 10
      END IF

      IF (mb(1)==mexpov .AND. ma(1)<mxexp1-ndig-2 .AND. mb(2)<0) THEN
        IF (krad==0) THEN
          CALL fmi2m(180,mc)
        ELSE
          CALL fmpi(mc)
        END IF
        IF (m01(2)<0) mc(2) = -mc(2)
        GO TO 10
      END IF

      IF (mb(1)==mexpun .AND. ma(2)==0) THEN
        IF (mb(2)<0) THEN
          IF (krad==0) THEN
            CALL fmi2m(180,mc)
          ELSE
            CALL fmpi(mc)
          END IF
        ELSE
          CALL fmi2m(0,mc)
        END IF
        GO TO 10
      END IF

      IF (mb(1)==mexpun .AND. (-ma(1))<mxexp1-ndig-2) THEN
        IF (krad==0) THEN
          CALL fmi2m(90,mc)
        ELSE
          CALL fmpi(mc)
          CALL fmdivi(mc,2,mc)
        END IF
        IF (m01(2)<0) mc(2) = -mc(2)
        GO TO 10
      END IF

!             Determine the quadrant for the result, then use FMATAN.

      IF (ma(2)>=0 .AND. mb(2)>0) jquad = 1
      IF (ma(2)>=0 .AND. mb(2)<0) jquad = 2
      IF (ma(2)<0 .AND. mb(2)<0) jquad = 3
      IF (ma(2)<0 .AND. mb(2)>0) jquad = 4

      CALL fmdiv(m01,m02,mc)
      mc(2) = abs(mc(2))
      CALL fmatan(mc,mc)

      IF (jquad==2 .OR. jquad==3) THEN
        IF (krad==0) THEN
          CALL fmi2m(180,m05)
          CALL fmsub(m05,mc,mc)
        ELSE
          CALL fmpi(m05)
          CALL fmsub(m05,mc,mc)
        END IF
      END IF

      IF ((jquad==3 .OR. jquad==4) .AND. mc(1)/=munkno) mc(2) = -mc(2)

!             Round the result and return.

10    IF (kflag==1) kflag = 0
      kwarn = kwrnsv
      macmax = nint((ndsave-1)*alogm2+log(real(abs(mc(2))+1))/0.69315)
      mc(0) = min(mc(0),macca,maccb,macmax)
      CALL fmexit(mc,mc,ndsave,mxsave,kasave,kovun)
      RETURN
    END SUBROUTINE fmatn2
    SUBROUTINE fmbig(ma)

!     MA = The biggest representable FM number using the current base
!          and precision.
!          The smallest positive number is then 1.0/MA.
!          Because of rounding, 1.0/(1.0/MA) will then overflow.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: j, n1
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmntr
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      namest(ncall) = 'FMBIG '

      IF (mblogs/=mbase) CALL fmcons
      kflag = 0
      n1 = ndig + 1
      DO 10 j = 2, n1
        ma(j) = mbase - 1
10    CONTINUE
      ma(1) = mxexp + 1
      ma(0) = nint(ndig*alogm2)

      IF (ntrace/=0) CALL fmntr(1,ma,ma,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmbig
    SUBROUTINE fmcat(ma,ncat)

!  NCAT is returned as the category of MA.  This is used by the various
!  arithmetic routines to handle special cases such as:
!  'number greater than 1' + 'underflowed result' is the first argument,
!  'overflowed result' / 'overflowed result' is 'unknown'.

!  NCAT       range

!   1.         -OV                OV stands for overflowed results.
!   2.   (-OV   , -OVTH)             ( MA(1) .GE. MAXEXP+2 )
!   3.   (-OVTH ,    -1)
!   4.         -1                 OVTH stands for a representable
!   5.   (-1    , -UNTH)               number near the overflow
!   6.   (-UNTH ,   -UN)               threshold.
!   7.         -UN                     ( MA(1) .GE. MAXEXP-NDIG+1 )
!   8.          0
!   9.         +UN                UN stands for underflowed results.
!  10.   (+UN   , +UNTH)             ( MA(1) .LE. -MAXEXP-1 )
!  11.   (+UNTH ,    +1)
!  12.         +1                 UNTH stands for a representable
!  13.   (+1    , +OVTH)               number near the underflow
!  14.   (+OVTH ,   +OV)               threshold.
!  15.         +OV                     ( MA(1) .LE. -MAXEXP+NDIG-1 )
!  16.       UNKNOWN

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, int
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ncat
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, mxexp1
      INTEGER :: j, nlast
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
!             Check for special symbols.

      ncat = 16
      IF (ma(1)==munkno) RETURN

      IF (ma(1)==mexpov) THEN
        ncat = 15
        IF (ma(2)<0) ncat = 1
        RETURN
      END IF

      IF (ma(1)==mexpun) THEN
        ncat = 9
        IF (ma(2)<0) ncat = 7
        RETURN
      END IF

      IF (ma(2)==0) THEN
        ncat = 8
        RETURN
      END IF

!             Check for +1 or -1.

      ma2 = abs(ma(2))
      IF (ma(1)==1 .AND. ma2==1) THEN
        nlast = ndig + 1
        IF (nlast>=3) THEN
          DO 10 j = 3, nlast
            IF (ma(j)/=0) GO TO 20
10        CONTINUE
        END IF
        ncat = 12
        IF (ma(2)<0) ncat = 4
        RETURN
      END IF

20    mxexp1 = int(mxexp2/2.01D0)
      IF (ma(1)>=mxexp1-ndig+1) THEN
        ncat = 14
        IF (ma(2)<0) ncat = 2
        RETURN
      END IF

      IF (ma(1)>=1) THEN
        ncat = 13
        IF (ma(2)<0) ncat = 3
        RETURN
      END IF

      IF (ma(1)>=-mxexp1+ndig) THEN
        ncat = 11
        IF (ma(2)<0) ncat = 5
        RETURN
      END IF

      IF (ma(1)>=-mxexp2) THEN
        ncat = 10
        IF (ma(2)<0) ncat = 6
        RETURN
      END IF

      RETURN
    END SUBROUTINE fmcat
    SUBROUTINE fmchsh(ma,mb,mc)

!  MB = COSH(MA),    MC = SINH(MA)

!  If both the hyperbolic sine and cosine are needed, this routine
!  is faster than calling both FMCOSH and FMSINH.

!  MB and MC must be distinct arrays.

      IMPLICIT NONE

!             Scratch array usage during FMCHSH:   M01 - M04

! .. Intrinsic Functions ..
      INTRINSIC abs, log, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, macca, macmax, mxsave
      INTEGER :: k, kasave, kovun, kreslt, kwrnsv, ncsave, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmcosh, fmentr, fmeq, fmeq2, fmexit, fmi2m, &
        fmntr, fmntrj, fmprnt, fmsinh, fmsqr, fmsqrt
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      macca = ma(0)
      ma2 = ma(2)
      IF (abs(ma(1))>mexpab) THEN
        ncsave = ncall
        CALL fmentr('FMCHSH',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
        IF (ma(1)==munkno) kovun = 2
        ncall = ncsave + 1
        CALL fmeq(ma,m04)
        m04(0) = nint(ndig*alogm2)
        m04(2) = abs(m04(2))
        CALL fmcosh(m04,mb)
        CALL fmsinh(m04,mc)
        GO TO 10
      ELSE
        ncall = ncall + 1
        namest(ncall) = 'FMCHSH'
        IF (ntrace/=0) CALL fmntr(2,ma,ma,1)
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
        ndsave = ndig
        IF (ncall==1) THEN
          k = max(ngrd52,2)
          ndig = max(ndig+k,2)
          IF (ndig>ndg2mx) THEN
            ncall = ncall - 1
            ndig = ndsave
            CALL fmeq(ma,m04)
            CALL fmcosh(m04,mb)
            CALL fmsinh(m04,mc)
            kflag = -9
            RETURN
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 0
        mxsave = mxexp
        mxexp = mxexp2
      END IF

      CALL fmeq2(ma,m04,ndsave,ndig,0)
      m04(0) = nint(ndig*alogm2)
      m04(2) = abs(m04(2))

      k = 1
      IF (m04(1)==0 .AND. m04(2)/=0) THEN
        IF (mbase/m04(2)>=100) k = 2
      END IF
      IF (m04(1)>=0 .AND. m04(2)/=0 .AND. k==1) THEN
        CALL fmcosh(m04,mb)
        IF (mb(1)>ndig) THEN
          CALL fmeq(mb,mc)
          GO TO 10
        END IF
        CALL fmsqr(mb,m03)
        CALL fmi2m(-1,m02)
        CALL fmadd(m03,m02,m03)
        CALL fmsqrt(m03,mc)
      ELSE
        CALL fmsinh(m04,mc)
        CALL fmsqr(mc,m03)
        CALL fmi2m(1,m02)
        CALL fmadd(m03,m02,m03)
        CALL fmsqrt(m03,mb)
      END IF

!             Round and return.

10    macmax = nint((ndsave-1)*alogm2+log(real(abs(mc(2))+1))/0.69315)
      mc(0) = min(mc(0),macca,macmax)
      IF (ma2<0 .AND. mc(1)/=munkno) mc(2) = -mc(2)
      CALL fmeq2(mc,mc,ndig,ndsave,1)
      macmax = nint((ndsave-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
      mb(0) = min(mb(0),macca,macmax)
      IF (kovun==2) THEN
        kwrnsv = kwarn
        kwarn = 0
      END IF
      CALL fmexit(mb,mb,ndsave,mxsave,kasave,kovun)
      IF (kovun==2) THEN
        kwarn = kwrnsv
      END IF
      IF (ntrace/=0) THEN
        IF (abs(ntrace)>=1 .AND. ncall+1<=lvltrc) THEN
          IF (ntrace<0) THEN
            CALL fmntrj(mc,ndig)
          ELSE
            CALL fmprnt(mc)
          END IF
        END IF
      END IF
      RETURN
    END SUBROUTINE fmchsh
    FUNCTION fmcomp(ma,lrel,mb)

!  Logical comparison of FM numbers MA and MB.

!  LREL is a CHARACTER *2 description of the comparison to be done:
!  LREL = 'EQ' returns FMCOMP = .TRUE. if MA.EQ.MB
!       = 'NE', 'GE', 'GT', 'LE', 'LT' also work like a logical IF.

!  For comparisons involving 'UNKNOWN' or two identical special symbols
!  such as +OVERFLOW,'EQ',+OVERFLOW, FMCOMP is returned FALSE and a
!  KFLAG = -4 error condition is returned.

!  Some compilers object to functions with side effects such as
!  changing KFLAG or other common variables.  Blocks of code that
!  modify common are identified by:
!      C                                                 DELETE START
!        ...
!      C                                                 DELETE STOP
!  These may be removed or commented out to produce a function without
!  side effects.  This disables trace printing in FMCOMP, and error
!  codes are not returned in KFLAG.

      IMPLICIT NONE

! .. Function Return Value ..
      LOGICAL :: fmcomp
! ..
! .. Intrinsic Functions ..
      INTRINSIC abs, int
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (2) :: lrel
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: j, jcomp, nlast
      CHARACTER (2) :: jrel
! ..
! .. External Subroutines ..
      EXTERNAL fmntrj, fmprnt
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
!                                                 DELETE START
      ncall = ncall + 1
      namest(ncall) = 'FMCOMP'

      IF (ncall<=lvltrc .AND. abs(ntrace)>=2) THEN
        WRITE (kw,90000)

        IF (ntrace>0) THEN
          CALL fmprnt(ma)
          WRITE (kw,90010) lrel
          CALL fmprnt(mb)
        ELSE
          CALL fmntrj(ma,ndig)
          WRITE (kw,90010) lrel
          CALL fmntrj(mb,ndig)
        END IF
      END IF
!                                                 DELETE STOP

!             JCOMP will be 1 if MA.GT.MB
!                           2 if MA.EQ.MB
!                           3 if MA.LT.MB

!             Check for special cases.

      jrel = lrel
      IF (lrel/='EQ' .AND. lrel/='NE' .AND. lrel/='LT' .AND. lrel/='GT' .AND. &
          lrel/='LE' .AND. lrel/='GE') THEN
        IF (lrel=='eq') THEN
          jrel = 'EQ'
        ELSE IF (lrel=='ne') THEN
          jrel = 'NE'
        ELSE IF (lrel=='lt') THEN
          jrel = 'LT'
        ELSE IF (lrel=='gt') THEN
          jrel = 'GT'
        ELSE IF (lrel=='le') THEN
          jrel = 'LE'
        ELSE IF (lrel=='ge') THEN
          jrel = 'GE'
        ELSE
          fmcomp = .FALSE.
!                                                 DELETE START
          kflag = -4
          IF (ncall/=1 .OR. kwarn<=0) GO TO 30
!                                                 DELETE STOP
          IF (kwarn<=0) GO TO 30
          WRITE (kw,90020) lrel
          IF (kwarn>=2) THEN
            STOP
          END IF
          GO TO 30
        END IF
      END IF

      IF (ma(1)==munkno .OR. mb(1)==munkno) THEN
        fmcomp = .FALSE.
!                                                 DELETE START
        kflag = -4
!                                                 DELETE STOP
        GO TO 30
      END IF

      IF (abs(ma(1))==mexpov .AND. ma(1)==mb(1) .AND. ma(2)==mb(2)) THEN
        fmcomp = .FALSE.
!                                                 DELETE START
        kflag = -4
        IF (ncall/=1 .OR. kwarn<=0) GO TO 30
!                                                 DELETE STOP
        IF (kwarn<=0) GO TO 30
        WRITE (kw,90030)
        IF (kwarn>=2) THEN
          STOP
        END IF
        GO TO 30
      END IF

!             Check for zero.

!                                                 DELETE START
      kflag = 0
!                                                 DELETE STOP
      IF (ma(2)==0) THEN
        jcomp = 2
        IF (mb(2)<0) jcomp = 1
        IF (mb(2)>0) jcomp = 3
        GO TO 20
      END IF
      IF (mb(2)==0) THEN
        jcomp = 1
        IF (ma(2)<0) jcomp = 3
        GO TO 20
      END IF
!             Check for opposite signs.

      IF (ma(2)>0 .AND. mb(2)<0) THEN
        jcomp = 1
        GO TO 20
      END IF
      IF (mb(2)>0 .AND. ma(2)<0) THEN
        jcomp = 3
        GO TO 20
      END IF

!             See which one is larger in absolute value.

      IF (ma(1)>mb(1)) THEN
        jcomp = 1
        GO TO 20
      END IF
      IF (mb(1)>ma(1)) THEN
        jcomp = 3
        GO TO 20
      END IF
      nlast = ndig + 1

      DO 10 j = 2, nlast
        IF (abs(ma(j))>abs(mb(j))) THEN
          jcomp = 1
          GO TO 20
        END IF
        IF (abs(mb(j))>abs(ma(j))) THEN
          jcomp = 3
          GO TO 20
        END IF
10    CONTINUE

      jcomp = 2

!             Now match the JCOMP value to the requested comparison.

20    IF (jcomp==1 .AND. ma(2)<0) THEN
        jcomp = 3
      ELSE IF (jcomp==3 .AND. mb(2)<0) THEN
        jcomp = 1
      END IF

      fmcomp = .FALSE.
      IF (jcomp==1 .AND. (jrel=='GT' .OR. jrel=='GE' .OR. jrel=='NE')) &
        fmcomp = .TRUE.

      IF (jcomp==2 .AND. (jrel=='EQ' .OR. jrel=='GE' .OR. jrel=='LE')) &
        fmcomp = .TRUE.

      IF (jcomp==3 .AND. (jrel=='NE' .OR. jrel=='LT' .OR. jrel=='LE')) &
        fmcomp = .TRUE.

30    CONTINUE
!                                                 DELETE START
      IF (ntrace/=0) THEN
        IF (ncall<=lvltrc .AND. abs(ntrace)>=1) THEN
          IF (kflag==0) THEN
            WRITE (kw,90040) ncall, int(mbase), ndig
          ELSE
            WRITE (kw,90050) ncall, int(mbase), ndig, kflag
          END IF
          IF (fmcomp) THEN
            WRITE (kw,90060)
          ELSE
            WRITE (kw,90070)
          END IF
        END IF
      END IF
      ncall = ncall - 1
!                                                 DELETE STOP
      RETURN
90000 FORMAT (' Input to FMCOMP')
90010 FORMAT (7X,'.',A2,'.')
90020 FORMAT (/' Error of type KFLAG = -4 in FM package in', &
        ' routine FMCOMP'//1X,A,' is not one of the six', &
        ' recognized comparisons.'//' .FALSE. has been',' returned.'/)
90030 FORMAT (/' Error of type KFLAG = -4 in FM package in routine', &
        ' FMCOMP'//' Two numbers in the same overflow or', &
        ' underflow category cannot be compared.'// &
        ' .FALSE. has been returned.'/)
90040 FORMAT (' FMCOMP',15X,'Call level =',I2,5X,'MBASE =',I10,5X,'NDIG =',I6)
90050 FORMAT (' FMCOMP',6X,'Call level =',I2,4X,'MBASE =',I10,4X,'NDIG =',I6, &
        4X,'KFLAG =',I3)
90060 FORMAT (7X,'.TRUE.')
90070 FORMAT (7X,'.FALSE.')
    END FUNCTION fmcomp
    SUBROUTINE fmcons

!  Set several saved machine precision constants.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC atan, dble, dint, int, log, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      mblogs = mbase
      alogmb = log(real(mbase))
      alogm2 = alogmb/log(2.0)
      alogmx = log(real(maxint))
      alogmt = alogmb/log(10.0)
      ngrd21 = int(2.0/alogmt+1.0)
      ngrd52 = int(5.0/alogmt+2.0)
      ngrd22 = int(2.0/alogmt+2.0)
      mexpab = dint(mxexp2/5)
      dlogmb = log(dble(mbase))
      dlogtn = log(10.0D0)
      dlogtw = log(2.0D0)
      dppi = 4.0D0*atan(1.0D0)
      dlogtp = log(2.0D0*dppi)
      dlogpi = log(dppi)
      dlogeb = -log(dpeps)/dlogmb

      RETURN
    END SUBROUTINE fmcons
    SUBROUTINE fmcos(ma,mb)

!  MB = COS(MA)

      IMPLICIT NONE

!             Scratch array usage during FMCOS:   M01 - M04

! .. Intrinsic Functions ..
      INTRINSIC abs, log, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: macca, macmax, mxsave
      INTEGER :: jcos, jsin, jswap, k, kasave, kovun, kreslt, ndsave, ndsv
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmcos2, fmdivi, fmentr, fmeq2, fmexit, fmi2m, fmmpy, &
        fmntr, fmpi, fmrdc, fmrslt, fmsin2, fmsqr, fmsqrt, fmsub, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (abs(ma(1))>mexpab .OR. ma(2)==0) THEN
        CALL fmentr('FMCOS ',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        namest(ncall) = 'FMCOS '
        IF (ntrace/=0) CALL fmntr(2,ma,ma,1)
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
        ndsave = ndig
        IF (ncall==1) THEN
          k = max(ngrd52,2)
          ndig = max(ndig+k,2)
          IF (ndig>ndg2mx) THEN
            kflag = -9
            CALL fmwarn
            ndig = ndsave
            kreslt = 12
            CALL fmrslt(ma,ma,mb,kreslt)
            IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
            ncall = ncall - 1
            RETURN
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 0
        mxsave = mxexp
        mxexp = mxexp2
      END IF

      macca = ma(0)
      CALL fmeq2(ma,mb,ndsave,ndig,0)
      mb(0) = nint(ndig*alogm2)
      mb(2) = abs(mb(2))

!             Reduce the argument, convert to radians if the input is
!             in degrees, and evaluate the function.

      CALL fmrdc(mb,mb,jsin,jcos,jswap)
      IF (mb(1)==munkno) GO TO 10
      IF (krad==0) THEN
        IF (mbspi/=mbase .OR. ndigpi<ndig) THEN
          ndsv = ndig
          ndig = min(ndig+2,ndg2mx)
          ncall = ncall + 1
          namest(ncall) = 'NOEQ  '
          CALL fmpi(mpisav)
          ncall = ncall - 1
          ndig = ndsv
        END IF
        CALL fmmpy(mb,mpisav,mb)
        CALL fmdivi(mb,180,mb)
      END IF
      IF (mb(1)/=munkno) THEN
        IF (jswap==0) THEN
          CALL fmcos2(mb,mb)
        ELSE
          IF (mb(1)<0 .OR. ndig<=50) THEN
            CALL fmsin2(mb,mb)
          ELSE
            CALL fmcos2(mb,mb)
            CALL fmi2m(1,m03)
            CALL fmsqr(mb,mb)
            CALL fmsub(m03,mb,mb)
            CALL fmsqrt(mb,mb)
          END IF
        END IF
      END IF

!             Append the sign, round, and return.

      IF (mb(1)/=munkno .AND. jcos==-1) mb(2) = -mb(2)
10    macmax = nint((ndsave-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
      mb(0) = min(mb(0),macca,macmax)
      CALL fmexit(mb,mb,ndsave,mxsave,kasave,kovun)
      RETURN
    END SUBROUTINE fmcos
    SUBROUTINE fmcos2(ma,mb)

!  Internal subroutine for MB = COS(MA) where 0.LE.MA.LE.1.

      IMPLICIT NONE

!             Scratch array usage during FMCOS2:   M01 - M04

!             LJSUMS = 8*(LUNPCK+1) allows for up to eight concurrent
!             sums.  Increasing this value will begin to improve the
!             speed of COS when the base is large and precision exceeds
!             about 1,500 decimal digits.

! .. Intrinsic Functions ..
      INTRINSIC int, log, max, min, nint, real, sqrt
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL :: alog2, alogt, b, t, tj
      REAL (KIND(0.0D0)) :: maxval
      INTEGER :: j, j2, k, k2, kpt, ktwo, kwrnsv, l, l2, large, n2, nbot, &
        ndsav1, ndsave, nterm
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmdivi, fmeq, fmeq2, fmi2m, fmipwr, fmmpy, &
        fmsqr, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mjsums(0:ljsums), &
        mlbsav(0:lunpck), mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), &
        mln4(0:lunpck), mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmsums/mjsums
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (ma(2)==0) THEN
        CALL fmi2m(1,mb)
        RETURN
      END IF
      ndsave = ndig
      kwrnsv = kwarn
      kwarn = 0

!             Use the direct series
!                  COS(X) = 1 - X**2/2! + X**4/4! - ...

!             The argument will be divided by 2**K2 before the series
!             is summed.  The series will be added as J2 concurrent
!             series.  The approximately optimal values of K2 and J2
!             are now computed to try to minimize the time required.
!             N2/2 is the approximate number of terms of the series
!             that will be needed, and L2 guard digits will be carried.

!             Since X is small when the series is summed, COS(X) - 1
!             is computed.  Then a version of the recovery formula can
!             be used that does not suffer from severe cancellation.

      b = real(mbase)
      k = ngrd52
      t = max(ndig-k,2)
      alog2 = log(2.0)
      alogt = log(t)
      tj = 0.03*alogmb*t**0.3333 + 1.85
      j2 = int(tj)
      j2 = max(1,min(j2,ljsums/ndg2mx))
      k2 = int(0.5*sqrt(t*alogmb/tj)+2.8)

      l = int(-(real(ma(1))*alogmb+log(real(ma(2))/b+ &
        real(ma(3))/(b*b)))/alog2-0.3)
      k2 = k2 - l
      IF (l<0) l = 0
      IF (k2<0) THEN
        k2 = 0
        j2 = int(.43*sqrt(t*alogmb/(alogt+real(l)*alog2))+.33)
      END IF
      IF (j2<=1) j2 = 1

      n2 = int(t*alogmb/(alogt+real(l)*alog2))
      l2 = int(log(real(n2)+2.0**k2)/alogmb)
      ndig = ndig + l2
      IF (ndig>ndg2mx) THEN
        kflag = -9
        CALL fmwarn
        mb(1) = munkno
        mb(2) = 1
        mb(0) = nint(ndig*alogm2)
        DO 10 j = 2, ndsave
          mb(j+1) = 0
10      CONTINUE
        ndig = ndsave
        kwarn = kwrnsv
        RETURN
      END IF
      ndsav1 = ndig

!             Divide the argument by 2**K2.

      CALL fmeq2(ma,m02,ndsave,ndig,0)
      ktwo = 1
      maxval = mxbase/2
      IF (k2>0) THEN
        DO 20 j = 1, k2
          ktwo = 2*ktwo
          IF (ktwo>maxval) THEN
            CALL fmdivi(m02,ktwo,m02)
            ktwo = 1
          END IF
20      CONTINUE
        IF (ktwo>1) CALL fmdivi(m02,ktwo,m02)
      END IF

!             Split into J2 concurrent sums and reduce NDIG while
!             computing each term in the sum as the terms get smaller.

      CALL fmsqr(m02,m02)
      CALL fmeq(m02,m03)
      m03(2) = -m03(2)
      nterm = 2
      DO 30 j = 1, j2
        nbot = nterm*(nterm-1)
        CALL fmdivi(m03,nbot,m03)
        nterm = nterm + 2
        kpt = (j-1)*(ndig+2)
        CALL fmeq(m03,mjsums(kpt))
        m03(2) = -m03(2)
30    CONTINUE
      IF (m02(1)<-ndig) GO TO 60
      CALL fmipwr(m02,j2,mb)

40    CALL fmmpy(m03,mb,m03)
      large = int(intmax/nterm)
      DO 50 j = 1, j2
        nbot = nterm*(nterm-1)
        IF (nterm>large .OR. nbot>mxbase) THEN
          CALL fmdivi(m03,nterm,m03)
          nbot = nterm - 1
          CALL fmdivi(m03,nbot,m03)
        ELSE
          CALL fmdivi(m03,nbot,m03)
        END IF
        kpt = (j-1)*(ndsav1+2)
        ndig = ndsav1
        CALL fmadd(mjsums(kpt),m03,mjsums(kpt))
        IF (kflag/=0) GO TO 60
        ndig = ndsav1 - int(mjsums(kpt+1)-m03(1))
        IF (ndig<2) ndig = 2
        m03(2) = -m03(2)
        nterm = nterm + 2
50    CONTINUE
      GO TO 40

!             Next put the J2 separate sums back together.

60    kflag = 0
      kpt = (j2-1)*(ndig+2)
      CALL fmeq(mjsums(kpt),mb)
      IF (j2>=2) THEN
        DO 70 j = 2, j2
          CALL fmmpy(m02,mb,mb)
          kpt = (j2-j)*(ndig+2)
          CALL fmadd(mb,mjsums(kpt),mb)
70      CONTINUE
      END IF

!             Reverse the effect of reducing the argument to
!             compute COS(MA).

      ndig = ndsav1
      IF (k2>0) THEN
        IF (ndsave<=20) THEN
          CALL fmi2m(2,m02)
          DO 80 j = 1, k2
            CALL fmadd(mb,m02,m03)
            CALL fmmpy(mb,m03,m03)
            CALL fmadd(m03,m03,mb)
80        CONTINUE
        ELSE
          DO 90 j = 1, k2
            CALL fmsqr(mb,m03)
            CALL fmadd(mb,mb,m02)
            CALL fmadd(m03,m02,m03)
            CALL fmadd(m03,m03,mb)
90        CONTINUE
        END IF
      END IF
      CALL fmi2m(1,m03)
      CALL fmadd(m03,mb,mb)

      CALL fmeq2(mb,mb,ndsav1,ndsave,1)
      ndig = ndsave
      kwarn = kwrnsv

      RETURN
    END SUBROUTINE fmcos2
    SUBROUTINE fmcosh(ma,mb)

!  MB = COSH(MA)

      IMPLICIT NONE

!             Scratch array usage during FMCOSH:   M01 - M03

! .. Intrinsic Functions ..
      INTRINSIC abs, int, log, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: macca, macmax, mxsave
      INTEGER :: k, kasave, kovun, kreslt, ndsave, nmethd
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmcsh2, fmdiv, fmdivi, fmentr, fmeq2, fmexit, &
        fmexp, fmi2m, fmntr, fmrslt, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (abs(ma(1))>mexpab) THEN
        CALL fmentr('FMCOSH',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        namest(ncall) = 'FMCOSH'
        IF (ntrace/=0) CALL fmntr(2,ma,ma,1)
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
        ndsave = ndig
        IF (ncall==1) THEN
          k = max(ngrd52,2)
          ndig = max(ndig+k,2)
          IF (ndig>ndg2mx) THEN
            kflag = -9
            CALL fmwarn
            ndig = ndsave
            kreslt = 12
            CALL fmrslt(ma,ma,mb,kreslt)
            IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
            ncall = ncall - 1
            RETURN
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 0
        mxsave = mxexp
        mxexp = mxexp2
      END IF

      macca = ma(0)
      CALL fmeq2(ma,mb,ndsave,ndig,0)
      mb(0) = nint(ndig*alogm2)
      mb(2) = abs(mb(2))
      IF (ma(2)==0) THEN
        CALL fmi2m(1,mb)
        GO TO 20
      END IF

!             Use a series for small arguments, FMEXP for large ones.

      IF (mb(1)==munkno) GO TO 20
      IF (mbase>99) THEN
        IF (mb(1)<=0) THEN
          nmethd = 1
        ELSE IF (mb(1)>=2) THEN
          nmethd = 2
        ELSE IF (abs(mb(2))<10) THEN
          nmethd = 1
        ELSE
          nmethd = 2
        END IF
      ELSE
        IF (mb(1)<=0) THEN
          nmethd = 1
        ELSE
          nmethd = 2
        END IF
      END IF

      IF (nmethd==2) GO TO 10
      CALL fmcsh2(mb,mb)
      GO TO 20

10    CALL fmexp(mb,mb)
      IF (mb(1)==mexpov) THEN
        GO TO 20
      ELSE IF (mb(1)==mexpun) THEN
        mb(1) = mexpov
        GO TO 20
      END IF
      IF (int(mb(1))<=(ndig+1)/2) THEN
        CALL fmi2m(1,m01)
        CALL fmdiv(m01,mb,m01)
        CALL fmadd(mb,m01,mb)
      END IF
      CALL fmdivi(mb,2,mb)

!             Round and return.

20    macmax = nint((ndsave-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
      mb(0) = min(mb(0),macca,macmax)
      CALL fmexit(mb,mb,ndsave,mxsave,kasave,kovun)
      RETURN
    END SUBROUTINE fmcosh
    SUBROUTINE fmcsh2(ma,mb)

!  Internal subroutine for MB = COSH(MA).

      IMPLICIT NONE

!             Scratch array usage during FMCSH2:   M01 - M03

!             LJSUMS = 8*(LUNPCK+1) allows for up to eight concurrent
!             sums.  Increasing this value will begin to improve the
!             speed of COSH when the base is large and precision exceeds
!             about 1,500 decimal digits.

! .. Intrinsic Functions ..
      INTRINSIC int, log, max, min, nint, real, sqrt
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL :: alog2, alogt, b, t, tj
      REAL (KIND(0.0D0)) :: maxval
      INTEGER :: j, j2, k, k2, kpt, ktwo, kwrnsv, l, l2, large, n2, nbot, &
        ndsav1, ndsave, nterm
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmdivi, fmeq, fmeq2, fmi2m, fmipwr, fmmpy, &
        fmsqr, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mjsums(0:ljsums), &
        mlbsav(0:lunpck), mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), &
        mln4(0:lunpck), mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmsums/mjsums
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (ma(2)==0) THEN
        CALL fmi2m(1,mb)
        RETURN
      END IF
      ndsave = ndig
      kwrnsv = kwarn
      kwarn = 0

!             Use the direct series
!                  COSH(X) = 1 + X**2/2! + X**4/4! - ...

!             The argument will be divided by 2**K2 before the series
!             is summed.  The series will be added as J2 concurrent
!             series.  The approximately optimal values of K2 and J2
!             are now computed to try to minimize the time required.
!             N2/2 is the approximate number of terms of the series
!             that will be needed, and L2 guard digits will be carried.

!             Since X is small when the series is summed, COSH(X) - 1
!             is computed.  Then a version of the recovery formula can
!             be used that does not suffer from severe cancellation.

      b = real(mbase)
      k = ngrd52
      t = max(ndig-k,2)
      alog2 = log(2.0)
      alogt = log(t)
      tj = 0.03*alogmb*t**0.3333 + 1.85
      j2 = int(tj)
      j2 = max(1,min(j2,ljsums/ndg2mx))
      k2 = int(0.5*sqrt(t*alogmb/tj)+2.8)

      l = int(-(real(ma(1))*alogmb+log(real(ma(2))/b+ &
        real(ma(3))/(b*b)))/alog2-0.3)
      k2 = k2 - l
      IF (l<0) l = 0
      IF (k2<0) THEN
        k2 = 0
        j2 = int(.43*sqrt(t*alogmb/(alogt+real(l)*alog2))+.33)
      END IF
      IF (j2<=1) j2 = 1

      n2 = int(t*alogmb/(alogt+real(l)*alog2))
      l2 = int(log(real(n2)+2.0**k2)/alogmb)
      ndig = ndig + l2
      IF (ndig>ndg2mx) THEN
        kflag = -9
        CALL fmwarn
        mb(1) = munkno
        mb(2) = 1
        mb(0) = nint(ndig*alogm2)
        DO 10 j = 2, ndsave
          mb(j+1) = 0
10      CONTINUE
        ndig = ndsave
        kwarn = kwrnsv
        RETURN
      END IF
      ndsav1 = ndig
      CALL fmeq2(ma,m02,ndsave,ndig,0)

!             Divide the argument by 2**K2.

      ktwo = 1
      maxval = mxbase/2
      IF (k2>0) THEN
        DO 20 j = 1, k2
          ktwo = 2*ktwo
          IF (ktwo>maxval) THEN
            CALL fmdivi(m02,ktwo,m02)
            ktwo = 1
          END IF
20      CONTINUE
        IF (ktwo>1) CALL fmdivi(m02,ktwo,m02)
      END IF

!             Split into J2 concurrent sums and reduce NDIG while
!             computing each term in the sum as the terms get smaller.

      CALL fmsqr(m02,m02)
      CALL fmeq(m02,m03)
      nterm = 2
      DO 30 j = 1, j2
        nbot = nterm*(nterm-1)
        CALL fmdivi(m03,nbot,m03)
        nterm = nterm + 2
        kpt = (j-1)*(ndig+2)
        CALL fmeq(m03,mjsums(kpt))
30    CONTINUE
      IF (m02(1)<-ndig) GO TO 60
      CALL fmipwr(m02,j2,mb)

40    CALL fmmpy(m03,mb,m03)
      large = int(intmax/nterm)
      DO 50 j = 1, j2
        nbot = nterm*(nterm-1)
        IF (nterm>large .OR. nbot>mxbase) THEN
          CALL fmdivi(m03,nterm,m03)
          nbot = nterm - 1
          CALL fmdivi(m03,nbot,m03)
        ELSE
          CALL fmdivi(m03,nbot,m03)
        END IF
        kpt = (j-1)*(ndsav1+2)
        ndig = ndsav1
        CALL fmadd(mjsums(kpt),m03,mjsums(kpt))
        IF (kflag/=0) GO TO 60
        ndig = ndsav1 - int(mjsums(kpt+1)-m03(1))
        IF (ndig<2) ndig = 2
        nterm = nterm + 2
50    CONTINUE
      GO TO 40

!             Next put the J2 separate sums back together.

60    kflag = 0
      kpt = (j2-1)*(ndig+2)
      CALL fmeq(mjsums(kpt),mb)
      IF (j2>=2) THEN
        DO 70 j = 2, j2
          CALL fmmpy(m02,mb,mb)
          kpt = (j2-j)*(ndig+2)
          CALL fmadd(mb,mjsums(kpt),mb)
70      CONTINUE
      END IF

!             Reverse the effect of reducing the argument to
!             compute COSH(MA).

      ndig = ndsav1
      IF (k2>0) THEN
        IF (ndsave<=20) THEN
          CALL fmi2m(2,m02)
          DO 80 j = 1, k2
            CALL fmadd(mb,m02,m03)
            CALL fmmpy(mb,m03,m03)
            CALL fmadd(m03,m03,mb)
80        CONTINUE
        ELSE
          DO 90 j = 1, k2
            CALL fmsqr(mb,m03)
            CALL fmadd(mb,mb,m02)
            CALL fmadd(m03,m02,m03)
            CALL fmadd(m03,m03,mb)
90        CONTINUE
        END IF
      END IF
      CALL fmi2m(1,m03)
      CALL fmadd(m03,mb,mb)

      CALL fmeq2(mb,mb,ndsav1,ndsave,1)
      ndig = ndsave
      kwarn = kwrnsv

      RETURN
    END SUBROUTINE fmcsh2
    SUBROUTINE fmcssn(ma,mb,mc)

!  MB = COS(MA),    MC = SIN(MA)

!  If both the sine and cosine are needed, this routine is faster
!  than calling both FMCOS and FMSIN.

!  MB and MC must be distinct arrays.

      IMPLICIT NONE

!             Scratch array usage during FMCSSN:   M01 - M05

! .. Intrinsic Functions ..
      INTRINSIC abs, log, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, macca, macmax, mxsave
      INTEGER :: jcos, jsin, jswap, k, kasave, kovun, kreslt, kwrnsv, ncsave, &
        ndsave, ndsv
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmcos, fmcos2, fmdivi, fmentr, fmeq, fmeq2, fmexit, &
        fmi2m, fmmpy, fmntr, fmntrj, fmpi, fmprnt, fmrdc, fmsin, fmsin2, &
        fmsqr, fmsqrt, fmsub
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      macca = ma(0)
      ma2 = ma(2)
      IF (abs(ma(1))>mexpab .OR. ma(2)==0) THEN
        ncsave = ncall
        CALL fmentr('FMCSSN',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
        IF (ma(1)==munkno) kovun = 2
        ncall = ncsave + 1
        CALL fmeq(ma,m05)
        m05(0) = nint(ndig*alogm2)
        m05(2) = abs(m05(2))
        CALL fmcos(m05,mb)
        CALL fmsin(m05,mc)
        GO TO 10
      ELSE
        ncall = ncall + 1
        namest(ncall) = 'FMCSSN'
        IF (ntrace/=0) CALL fmntr(2,ma,ma,1)
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
        ndsave = ndig
        IF (ncall==1) THEN
          k = max(ngrd52,2)
          ndig = max(ndig+k,2)
          IF (ndig>ndg2mx) THEN
            ncall = ncall - 1
            ndig = ndsave
            CALL fmeq(ma,m05)
            CALL fmcos(m05,mb)
            CALL fmsin(m05,mc)
            kflag = -9
            RETURN
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 0
        mxsave = mxexp
        mxexp = mxexp2
      END IF

      IF (ma(2)==0) THEN
        CALL fmi2m(1,mb)
        CALL fmi2m(0,mc)
        GO TO 10
      END IF

      CALL fmeq2(ma,mb,ndsave,ndig,0)
      mb(0) = nint(ndig*alogm2)
      mb(2) = abs(mb(2))

!             Reduce the argument, convert to radians if the input is
!             in degrees, and evaluate the functions.

      CALL fmrdc(mb,mb,jsin,jcos,jswap)
      IF (mb(1)==munkno) THEN
        CALL fmeq(mb,mc)
        GO TO 10
      END IF
      IF (krad==0) THEN
        IF (mbspi/=mbase .OR. ndigpi<ndig) THEN
          ndsv = ndig
          ndig = min(ndig+2,ndg2mx)
          ncall = ncall + 1
          namest(ncall) = 'NOEQ  '
          CALL fmpi(mpisav)
          ncall = ncall - 1
          ndig = ndsv
        END IF
        CALL fmmpy(mb,mpisav,mb)
        CALL fmdivi(mb,180,mb)
      END IF
      IF (mb(1)/=munkno) THEN
        IF (jswap==0) THEN
          IF (mb(1)<0) THEN
            CALL fmsin2(mb,mc)
            mc(2) = jsin*mc(2)
            CALL fmsqr(mc,m03)
            CALL fmi2m(1,m02)
            CALL fmsub(m02,m03,m03)
            CALL fmsqrt(m03,mb)
            mb(2) = jcos*mb(2)
          ELSE
            CALL fmcos2(mb,mb)
            mb(2) = jcos*mb(2)
            CALL fmsqr(mb,m03)
            CALL fmi2m(1,m02)
            CALL fmsub(m02,m03,m03)
            CALL fmsqrt(m03,mc)
            mc(2) = jsin*mc(2)
          END IF
        ELSE
          IF (mb(1)<0) THEN
            CALL fmsin2(mb,mb)
            mb(2) = jcos*mb(2)
            CALL fmsqr(mb,m03)
            CALL fmi2m(1,m02)
            CALL fmsub(m02,m03,m03)
            CALL fmsqrt(m03,mc)
            mc(2) = jsin*mc(2)
          ELSE
            CALL fmcos2(mb,mc)
            mc(2) = jsin*mc(2)
            CALL fmsqr(mc,m03)
            CALL fmi2m(1,m02)
            CALL fmsub(m02,m03,m03)
            CALL fmsqrt(m03,mb)
            mb(2) = jcos*mb(2)
          END IF
        END IF
      ELSE
        CALL fmeq(mb,mc)
      END IF

!             Round and return.

10    macmax = nint((ndsave-1)*alogm2+log(real(abs(mc(2))+1))/0.69315)
      mc(0) = min(mc(0),macca,macmax)
      IF (ma2<0 .AND. mc(1)/=munkno) mc(2) = -mc(2)
      CALL fmeq2(mc,mc,ndig,ndsave,1)
      mb(0) = min(mb(0),macca,macmax)
      IF (kovun==2) THEN
        kwrnsv = kwarn
        kwarn = 0
      END IF
      CALL fmexit(mb,mb,ndsave,mxsave,kasave,kovun)
      IF (kovun==2) THEN
        kwarn = kwrnsv
      END IF
      IF (ntrace/=0) THEN
        IF (abs(ntrace)>=1 .AND. ncall+1<=lvltrc) THEN
          IF (ntrace<0) THEN
            CALL fmntrj(mc,ndig)
          ELSE
            CALL fmprnt(mc)
          END IF
        END IF
      END IF
      RETURN
    END SUBROUTINE fmcssn
    SUBROUTINE fmdbl(a,b,c)

!  C = A + B.  All are double precision.  This routine tries to
!  force the compiler to round C to double precision accuracy.
!  Some compilers allow double precision loops like the ones in
!  FMSET and FMDM to be done in extended precision, which defeats
!  the routine's attempt to determine double precision accuracy.
!  This can lead to doing too few Newton steps and failing to
!  get sufficient accuracy in several FM routines.

! .. Scalar Arguments ..
      REAL (KIND(0.0D0)) :: a, b, c
! ..
      c = a + b
      RETURN
    END SUBROUTINE fmdbl
    SUBROUTINE fmdig(nstack,kst)

!  Compute the number of intermediate digits to be used in Newton
!  iteration.  This assumes that a starting approximation that is
!  accurate to double precision is used, and the root is simple.

!  KST is the number of iterations needed for final accuracy NDIG.
!  NSTACK(J) holds the value of NDIG to be used for the
!            Jth iteration.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC int
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kst
! ..
! .. Array Arguments ..
      INTEGER :: nstack(19)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: y
      INTEGER :: j, jt, l, nd, ndt, ne
! ..
! .. External Subroutines ..
      EXTERNAL fmcons
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons

!             NE is the maximum number of base MBASE digits that
!             can be used in the first Newton iteration.

      ne = int(1.9D0*dlogeb)

!             Fill the intermediate digit stack (backwards).

      kst = 1
      nd = ndig
      nstack(1) = nd
      IF (nd<ne .OR. nd<=2) RETURN

10    y = nd

!             The 1.9 accounts for the fact that the number of correct
!             digits approximately doubles at each iteration.

      ndt = int(y/1.9D0)
      IF (2*ndt<=nd) ndt = ndt + 1
      nd = ndt
      kst = kst + 1
      nstack(kst) = nd
      IF (nd>ne .AND. nd>2) GO TO 10

!             Reverse the stack.

      l = kst/2
      DO 20 j = 1, l
        jt = nstack(j)
        nstack(j) = nstack(kst+1-j)
        nstack(kst+1-j) = jt
20    CONTINUE

      RETURN
    END SUBROUTINE fmdig
    SUBROUTINE fmdim(ma,mb,mc)

!  MC = DIM(MA,MB)

!  Positive difference.  MC = MA - MB  if MA.GE.MB,
!                           = 0        otherwise.

      IMPLICIT NONE

!             Scratch array usage during FMDIM:   M01 - M02

! .. Intrinsic Functions ..
      INTRINSIC abs, log, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: macca, maccb, macmax, mxsave
      INTEGER :: k, kasave, kovun, kreslt, kwrnsv, ndsave
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmentr, fmeq2, fmexit, fmi2m, fmntr, fmrslt, fmsub, &
        fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (abs(ma(1))>mexpab .OR. abs(mb(1))>mexpab) THEN
        CALL fmentr('FMDIM ',ma,mb,2,mc,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        namest(ncall) = 'FMDIM '
        IF (ntrace/=0) CALL fmntr(2,ma,mb,2)
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
        IF (mb(1)==mexpov .OR. mb(1)==mexpun) kovun = 1
        ndsave = ndig
        IF (ncall==1) THEN
          k = max(ngrd52-1,2)
          ndig = max(ndig+k,2)
          IF (ndig>ndg2mx) THEN
            kflag = -9
            CALL fmwarn
            ndig = ndsave
            kreslt = 12
            CALL fmrslt(ma,mb,mc,kreslt)
            IF (ntrace/=0) CALL fmntr(1,mc,mc,1)
            ncall = ncall - 1
            RETURN
          END IF
        END IF
        kasave = kaccsw
        mxsave = mxexp
        mxexp = mxexp2
      END IF
      kwrnsv = kwarn
      kwarn = 0
      mxexp = mxsave

      macca = ma(0)
      maccb = mb(0)
      CALL fmeq2(ma,m01,ndsave,ndig,0)
      m01(0) = nint(ndig*alogm2)
      CALL fmeq2(mb,m02,ndsave,ndig,0)
      m02(0) = nint(ndig*alogm2)

      IF (fmcomp(m01,'LT',m02)) THEN
        CALL fmi2m(0,mc)
      ELSE
        CALL fmsub(m01,m02,mc)
      END IF

      IF (kflag==1) kflag = 0
      kwarn = kwrnsv
      macmax = nint((ndsave-1)*alogm2+log(real(abs(mc(2))+1))/0.69315)
      mc(0) = min(mc(0),macca,maccb,macmax)
      CALL fmexit(mc,mc,ndsave,mxsave,kasave,kovun)
      RETURN
    END SUBROUTINE fmdim
    SUBROUTINE fmdiv(ma,mb,mc)

!  MC = MA / MB

!  This routine performs the trace printing for division.
!  FMDIV2 is used to do the arithmetic.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. External Subroutines ..
      EXTERNAL fmdiv2, fmntr
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (ntrace/=0) THEN
        namest(ncall) = 'FMDIV '
        CALL fmntr(2,ma,mb,2)

        CALL fmdiv2(ma,mb,mc)

        CALL fmntr(1,mc,mc,1)
      ELSE
        CALL fmdiv2(ma,mb,mc)
      END IF
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmdiv
    SUBROUTINE fmdiv2(ma,mb,mc)

!  Internal division routine.  MC = MA / MB

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dble, dint, int, log, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, ma2p, macca, maccb, maxmwa, mb1, mb2, mb2p, mbm1, &
        mcarry, md2b, mkt, mlmax, mlr, mqd
      REAL (KIND(0.0D0)) :: xb, xbase, xbr, xmwa
      INTEGER :: j, jb, jl, ka, kb, kl, kptmwa, kreslt, n1, ng, nguard, nl, &
        nmbwds, nzdmb
! ..
! .. External Subroutines ..
      EXTERNAL fmargs, fmcons, fmim, fmmove, fmrnd, fmrslt, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      macca = ma(0)
      maccb = mb(0)
      IF (abs(ma(1))>mexpab .OR. abs(mb(1))>mexpab .OR. kdebug==1) THEN
        CALL fmargs('FMDIV ',2,ma,mb,kreslt)
        IF (kreslt/=0) THEN
          ncall = ncall + 1
          namest(ncall) = 'FMDIV '
          CALL fmrslt(ma,mb,mc,kreslt)
          ncall = ncall - 1
          RETURN
        END IF
      ELSE
        IF (mb(2)==0) THEN
          CALL fmim(0,mc)
          mc(1) = munkno
          mc(2) = 1
          mc(0) = nint(ndig*alogm2)
          namest(ncall) = 'FMDIV '
          kflag = -4
          CALL fmwarn
          RETURN
        END IF
        IF (ma(2)==0) THEN
          CALL fmim(0,mc)
          mc(0) = min(macca,maccb)
          RETURN
        END IF
      END IF
      kflag = 0

!             NGUARD is the number of guard digits used.

      IF (ncall>1) THEN
        nguard = ngrd21
        IF (nguard>ndig) nguard = ndig
      ELSE
        nguard = ngrd52 - 1
      END IF
      ma2p = abs(ma(2))
      mb2p = abs(mb(2))
      n1 = ndig + 1
      ng = ndig + nguard

!             Copy MA into the working array.

      DO 10 j = 3, n1
        mwa(j+1) = ma(j)
10    CONTINUE
      mwa(1) = ma(1) - mb(1) + 1
      mwa(2) = 0
      nl = n1 + nguard + 3
      DO 20 j = ndig + 3, nl
        mwa(j) = 0
20    CONTINUE

!             Save the sign of MA and MB and then work only with
!             positive numbers.

      ma2 = ma(2)
      mb1 = mb(1)
      mb2 = mb(2)
      ma(2) = ma2p
      mwa(3) = ma(2)
      mb(1) = 0
      mb(2) = mb2p

!             NMBWDS is the number of words of MB used to
!             compute the estimated quotient digit MQD.

      nmbwds = 4
      IF (mbase<100) nmbwds = 7

!             XB is an approximation of MB used in
!             estimating the quotient digits.

      xbase = dble(mbase)
      xb = 0
      jl = nmbwds
      IF (jl<=n1) THEN
        DO 30 j = 2, jl
          xb = xb*xbase + dble(mb(j))
30      CONTINUE
      ELSE
        DO 40 j = 2, jl
          IF (j<=n1) THEN
            xb = xb*xbase + dble(mb(j))
          ELSE
            xb = xb*xbase
          END IF
40      CONTINUE
      END IF
      IF (jl+1<=n1) xb = xb + dble(mb(jl+1))/xbase
      xbr = 1.0D0/xb

!             MLMAX determines when to normalize all of MWA.

      mbm1 = mbase - 1
      mlmax = maxint/mbm1
      mkt = intmax - mbase
      mlmax = min(mlmax,mkt)

!             Count the trailing zero digits of MB.

      DO 50 j = n1, 2, -1
        IF (mb(j)/=0) THEN
          nzdmb = n1 - j
          GO TO 60
        END IF
50    CONTINUE

!             MAXMWA is an upper bound on the size of values in MWA
!             divided by MBASE-1.  It is used to determine whether
!             normalization can be postponed.

60    maxmwa = 0

!             KPTMWA points to the next digit in the quotient.

      kptmwa = 2

!             This is the start of the division loop.

!             XMWA is an approximation of the active part of MWA
!             used in estimating quotient digits.

70    kl = kptmwa + nmbwds - 1
      IF (kl<=nl) THEN
        xmwa = ((dble(mwa(kptmwa))*xbase+dble(mwa(kptmwa+1)))*xbase+dble(mwa( &
          kptmwa+2)))*xbase + dble(mwa(kptmwa+3))
        DO 80 j = kptmwa + 4, kl
          xmwa = xmwa*xbase + dble(mwa(j))
80      CONTINUE
      ELSE
        xmwa = dble(mwa(kptmwa))
        DO 90 j = kptmwa + 1, kl
          IF (j<=nl) THEN
            xmwa = xmwa*xbase + dble(mwa(j))
          ELSE
            xmwa = xmwa*xbase
          END IF
90      CONTINUE
      END IF

!             MQD is the estimated quotient digit.

      mqd = dint(xmwa*xbr)
      IF (mqd<0) mqd = mqd - 1

      IF (mqd>0) THEN
        maxmwa = maxmwa + mqd
      ELSE
        maxmwa = maxmwa - mqd
      END IF

!             See if MWA must be normalized.

      ka = kptmwa + 1
      kb = min(ka+ndig-1-nzdmb,nl)
      IF (maxmwa>=mlmax) THEN
        DO 100 j = kb, ka, -1
          IF (mwa(j)<0) THEN
            mcarry = int((-mwa(j)-1)/mbase) + 1
            mwa(j) = mwa(j) + mcarry*mbase
            mwa(j-1) = mwa(j-1) - mcarry
          ELSE IF (mwa(j)>=mbase) THEN
            mcarry = -int(mwa(j)/mbase)
            mwa(j) = mwa(j) + mcarry*mbase
            mwa(j-1) = mwa(j-1) - mcarry
          END IF
100     CONTINUE
        xmwa = 0
        IF (kl<=nl) THEN
          DO 110 j = kptmwa, kl
            xmwa = xmwa*xbase + dble(mwa(j))
110       CONTINUE
        ELSE
          DO 120 j = kptmwa, kl
            IF (j<=nl) THEN
              xmwa = xmwa*xbase + dble(mwa(j))
            ELSE
              xmwa = xmwa*xbase
            END IF
120       CONTINUE
        END IF
        mqd = dint(xmwa*xbr)
        IF (mqd<0) mqd = mqd - 1
        IF (mqd>0) THEN
          maxmwa = mqd
        ELSE
          maxmwa = -mqd
        END IF
      END IF

!             Subtract MQD*MB from MWA.

      jb = ka - 2
      IF (mqd/=0) THEN

!             Major (Inner Loop)

        DO 130 j = ka, kb
          mwa(j) = mwa(j) - mqd*mb(j-jb)
130     CONTINUE
      END IF

      mwa(ka) = mwa(ka) + mwa(ka-1)*mbase
      mwa(kptmwa) = mqd

      kptmwa = kptmwa + 1
      IF (kptmwa<=ng) GO TO 70
      IF (mwa(2)==0 .AND. kptmwa<=ng+1) GO TO 70

      kl = kptmwa + nmbwds - 1
      IF (kl<=nl) THEN
        xmwa = ((dble(mwa(kptmwa))*xbase+dble(mwa(kptmwa+1)))*xbase+dble(mwa( &
          kptmwa+2)))*xbase + dble(mwa(kptmwa+3))
        DO 140 j = kptmwa + 4, kl
          xmwa = xmwa*xbase + dble(mwa(j))
140     CONTINUE
      ELSE
        xmwa = dble(mwa(kptmwa))
        DO 150 j = kptmwa + 1, kl
          IF (j<=nl) THEN
            xmwa = xmwa*xbase + dble(mwa(j))
          ELSE
            xmwa = xmwa*xbase
          END IF
150     CONTINUE
      END IF
      mqd = dint(xmwa*xbr)
      IF (mqd<0) mqd = mqd - 1
      mwa(kptmwa) = mqd
      mwa(kptmwa+1) = 0
      mwa(kptmwa+2) = 0

!             Final normalization.

      DO 160 j = kptmwa, 3, -1
        IF (mwa(j)<0) THEN
          mcarry = int((-mwa(j)-1)/mbase) + 1
          mwa(j) = mwa(j) + mcarry*mbase
          mwa(j-1) = mwa(j-1) - mcarry
        ELSE IF (mwa(j)>=mbase) THEN
          mcarry = -int(mwa(j)/mbase)
          mwa(j) = mwa(j) + mcarry*mbase
          mwa(j-1) = mwa(j-1) - mcarry
        END IF
160   CONTINUE

!             Round, affix the sign, and return.

      ma(2) = ma2
      mb(1) = mb1
      mb(2) = mb2
      IF (mwa(2)==0) THEN
        mlr = 2*mwa(ndig+3) + 1
        IF (mlr>=mbase) THEN
          IF (mlr-1>mbase .AND. mwa(n1+1)<mbase-1) THEN
            IF (kround/=0 .OR. ncall>1) THEN
              mwa(n1+1) = mwa(n1+1) + 1
              mwa(n1+2) = 0
            END IF
          ELSE
            CALL fmrnd(mwa,ndig,nguard,1)
          END IF
        END IF
      ELSE
        mlr = 2*mwa(ndig+2) + 1
        IF (mlr>=mbase) THEN
          IF (mlr-1>mbase .AND. mwa(n1)<mbase-1) THEN
            IF (kround/=0 .OR. ncall>1) THEN
              mwa(n1) = mwa(n1) + 1
              mwa(n1+1) = 0
            END IF
          ELSE
            CALL fmrnd(mwa,ndig,nguard,0)
          END IF
        END IF
      END IF
      CALL fmmove(mwa,mc)

      IF (kflag<0) THEN
        namest(ncall) = 'FMDIV '
        CALL fmwarn
      END IF

      IF (ma2*mb2<0) mc(2) = -mc(2)

      IF (kaccsw==1) THEN
        md2b = nint((ndig-1)*alogm2+log(real(abs(mc(2))+1))/0.69315)
        mc(0) = min(macca,maccb,md2b)
      ELSE
        mc(0) = min(macca,maccb)
      END IF
      RETURN
    END SUBROUTINE fmdiv2
    SUBROUTINE fmdivd(ma,mb,mc,md,me)

!  Double division routine.  MD = MA / MC,   ME = MB / MC

!  It is usually slightly faster to do two divisions that
!  have a common denominator with one call.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dble, dint, int, log, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck), &
        me(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, ma2p, macca, maccb, maccc, maxmwa, mb2, mb2p, mbm1, &
        mc1, mc2, mc2p, mcarry, md2b, mkt, mlmax, mlr, mqdmwa, mqdmwd, mtemp
      REAL (KIND(0.0D0)) :: xb, xbase, xbr, xmwa, xmwd
      INTEGER :: j, jb, jl, ka, kb, kl, kovun, kptmw, n1, ng, nguard, nl, &
        nmbwds, nzdmb
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmdiv2, fmeq, fmim, fmmove, fmntr, fmntrj, fmprnt, &
        fmrnd, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa), mwd(lmwa), mwe(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
      COMMON /fmwa/mwd, mwe
! ..
      ncall = ncall + 1
      IF (ntrace/=0) THEN
        namest(ncall) = 'FMDIVD'
        CALL fmntr(2,ma,mb,2)
        IF (abs(ntrace)>=2 .AND. ncall<=lvltrc) THEN
          IF (ntrace<0) THEN
            CALL fmntrj(mc,ndig)
          ELSE
            CALL fmprnt(mc)
          END IF
        END IF
      END IF

      IF (mblogs/=mbase) CALL fmcons
      macca = ma(0)
      maccb = mb(0)
      maccc = mc(0)
      IF (abs(ma(1))>mexpab .OR. abs(mb(1))>mexpab .OR. abs(mc(1))>mexpab) &
          THEN
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun .OR. mb(1)==mexpov .OR. &
          mb(1)==mexpun .OR. mc(1)==mexpov .OR. mc(1)==mexpun) kovun = 1
        IF (ma(1)==munkno .OR. mb(1)==munkno .OR. mc(1)==munkno) kovun = 2
        ncall = ncall + 1
        CALL fmdiv2(ma,mc,mwd)
        kb = kflag
        CALL fmdiv2(mb,mc,me)
        ncall = ncall - 1
        IF (((kflag<0 .OR. kb<0) .AND. kovun==0) .OR. ((kflag==-4 .OR. kb== &
            -4) .AND. kovun==1)) THEN
          IF (kflag==-4 .OR. kb==-4) THEN
            kflag = -4
          ELSE IF (kflag==-5 .OR. kb==-5) THEN
            kflag = -5
          ELSE
            kflag = min(kflag,kb)
          END IF
          namest(ncall) = 'FMDIVD'
          CALL fmwarn
        END IF
        CALL fmeq(mwd,md)
        GO TO 170
      END IF
      IF (mc(2)==0) THEN
        CALL fmim(0,md)
        md(1) = munkno
        md(2) = 1
        md(0) = nint(ndig*alogm2)
        CALL fmim(0,me)
        me(1) = munkno
        me(2) = 1
        me(0) = nint(ndig*alogm2)
        namest(ncall) = 'FMDIVD'
        kflag = -4
        CALL fmwarn
        GO TO 170
      END IF
      IF (ma(2)==0 .OR. mb(2)==0) THEN
        CALL fmdiv2(ma,mc,mwd)
        CALL fmdiv2(mb,mc,me)
        CALL fmeq(mwd,md)
        GO TO 170
      END IF
      kflag = 0

!             NGUARD is the number of guard digits used.

      IF (ncall>1) THEN
        nguard = ngrd21
        IF (nguard>ndig) nguard = ndig
      ELSE
        nguard = ngrd52 - 1
      END IF
      ma2p = abs(ma(2))
      mb2p = abs(mb(2))
      mc2p = abs(mc(2))
      IF ((mc2p>=ma2p .OR. mc2p>=mb2p) .AND. nguard<2) nguard = 2
      n1 = ndig + 1
      ng = ndig + nguard

!             Copy MA and MB into the working arrays.

      DO 10 j = 3, n1
        mwa(j+1) = ma(j)
        mwd(j+1) = mb(j)
10    CONTINUE
      mwa(1) = ma(1) - mc(1) + 1
      mwd(1) = mb(1) - mc(1) + 1
      mwa(2) = 0
      mwd(2) = 0
      nl = n1 + nguard + 3
      DO 20 j = ndig + 3, nl
        mwa(j) = 0
        mwd(j) = 0
20    CONTINUE

!             Save the signs and then work only with
!             positive numbers.

      ma2 = ma(2)
      mb2 = mb(2)
      mc1 = mc(1)
      mc2 = mc(2)
      ma(2) = ma2p
      mb(2) = mb2p
      mwa(3) = ma(2)
      mwd(3) = mb(2)
      mc(1) = 0
      mc(2) = mc2p

!             NMBWDS is the number of words used to compute
!             the estimated quotient digits.

      nmbwds = 4
      IF (mbase<100) nmbwds = 7

!             XB is an approximation of MC used in selecting
!             estimated quotients.

      xbase = dble(mbase)
      xb = 0
      jl = nmbwds
      IF (jl<=n1) THEN
        DO 30 j = 2, jl
          xb = xb*xbase + dble(mc(j))
30      CONTINUE
      ELSE
        DO 40 j = 2, jl
          IF (j<=n1) THEN
            xb = xb*xbase + dble(mc(j))
          ELSE
            xb = xb*xbase
          END IF
40      CONTINUE
      END IF
      IF (jl+1<=n1) xb = xb + dble(mc(jl+1))/xbase
      xbr = 1.0D0/xb

!             MLMAX determines when to normalize all of MWA.

      mbm1 = mbase - 1
      mlmax = maxint/mbm1
      mkt = intmax - mbase
      mlmax = min(mlmax,mkt)

!             Count the trailing zero digits of MC.

      DO 50 j = n1, 2, -1
        IF (mc(j)/=0) THEN
          nzdmb = n1 - j
          GO TO 60
        END IF
50    CONTINUE

!             MAXMWA is an upper bound on the size of values in MWA
!             divided by MBASE-1.  It is used to determine whether
!             normalization can be postponed.

60    maxmwa = 0

!             KPTMW points to the next digit in the quotient.

      kptmw = 2

!             This is the start of the division loop.

!             XMWA is an approximation of the active part of MWA
!             used in selecting estimated quotients.

70    kl = kptmw + nmbwds - 1
      IF (kl<=nl) THEN
        xmwa = ((dble(mwa(kptmw))*xbase+dble(mwa(kptmw+1)))*xbase+dble(mwa( &
          kptmw+2)))*xbase + dble(mwa(kptmw+3))
        xmwd = ((dble(mwd(kptmw))*xbase+dble(mwd(kptmw+1)))*xbase+dble(mwd( &
          kptmw+2)))*xbase + dble(mwd(kptmw+3))
        DO 80 j = kptmw + 4, kl
          xmwa = xmwa*xbase + dble(mwa(j))
          xmwd = xmwd*xbase + dble(mwd(j))
80      CONTINUE
      ELSE
        xmwa = dble(mwa(kptmw))
        xmwd = dble(mwd(kptmw))
        DO 90 j = kptmw + 1, kl
          IF (j<=nl) THEN
            xmwa = xmwa*xbase + dble(mwa(j))
            xmwd = xmwd*xbase + dble(mwd(j))
          ELSE
            xmwa = xmwa*xbase
            xmwd = xmwd*xbase
          END IF
90      CONTINUE
      END IF

!             MQDMWA and MQDMWD are the estimated quotient digits.

      mqdmwa = dint(xmwa*xbr)
      IF (mqdmwa<0) mqdmwa = mqdmwa - 1
      mqdmwd = dint(xmwd*xbr)
      IF (mqdmwd<0) mqdmwd = mqdmwd - 1

      maxmwa = maxmwa + max(abs(mqdmwa),abs(mqdmwd))

!             See if MWA and MWD must be normalized.

      ka = kptmw + 1
      kb = min(ka+ndig-1-nzdmb,nl)
      IF (maxmwa>=mlmax) THEN
        DO 100 j = kb, ka, -1
          IF (mwa(j)<0) THEN
            mcarry = int((-mwa(j)-1)/mbase) + 1
            mwa(j) = mwa(j) + mcarry*mbase
            mwa(j-1) = mwa(j-1) - mcarry
          ELSE IF (mwa(j)>=mbase) THEN
            mcarry = -int(mwa(j)/mbase)
            mwa(j) = mwa(j) + mcarry*mbase
            mwa(j-1) = mwa(j-1) - mcarry
          END IF
          IF (mwd(j)<0) THEN
            mcarry = int((-mwd(j)-1)/mbase) + 1
            mwd(j) = mwd(j) + mcarry*mbase
            mwd(j-1) = mwd(j-1) - mcarry
          ELSE IF (mwd(j)>=mbase) THEN
            mcarry = -int(mwd(j)/mbase)
            mwd(j) = mwd(j) + mcarry*mbase
            mwd(j-1) = mwd(j-1) - mcarry
          END IF
100     CONTINUE
        xmwa = 0
        xmwd = 0
        IF (kl<=nl) THEN
          DO 110 j = kptmw, kl
            xmwa = xmwa*xbase + dble(mwa(j))
            xmwd = xmwd*xbase + dble(mwd(j))
110       CONTINUE
        ELSE
          DO 120 j = kptmw, kl
            IF (j<=nl) THEN
              xmwa = xmwa*xbase + dble(mwa(j))
              xmwd = xmwd*xbase + dble(mwd(j))
            ELSE
              xmwa = xmwa*xbase
              xmwd = xmwd*xbase
            END IF
120       CONTINUE
        END IF
        mqdmwa = dint(xmwa*xbr)
        IF (mqdmwa<0) mqdmwa = mqdmwa - 1
        mqdmwd = dint(xmwd*xbr)
        IF (mqdmwd<0) mqdmwd = mqdmwd - 1
        maxmwa = max(abs(mqdmwa),abs(mqdmwd))
      END IF

!             Subtract MQDMWA*MC from MWA and MQDMWD*MC from MWD.

      jb = ka - 2

!             Major (Inner Loop)

      DO 130 j = ka, kb
        mtemp = mc(j-jb)
        mwa(j) = mwa(j) - mqdmwa*mtemp
        mwd(j) = mwd(j) - mqdmwd*mtemp
130   CONTINUE

      mwa(ka) = mwa(ka) + mwa(ka-1)*mbase
      mwd(ka) = mwd(ka) + mwd(ka-1)*mbase
      mwa(kptmw) = mqdmwa
      mwd(kptmw) = mqdmwd

      kptmw = kptmw + 1
      IF (kptmw<=ng) GO TO 70

      kl = kptmw + nmbwds - 1
      IF (kl<=nl) THEN
        xmwa = ((dble(mwa(kptmw))*xbase+dble(mwa(kptmw+1)))*xbase+dble(mwa( &
          kptmw+2)))*xbase + dble(mwa(kptmw+3))
        xmwd = ((dble(mwd(kptmw))*xbase+dble(mwd(kptmw+1)))*xbase+dble(mwd( &
          kptmw+2)))*xbase + dble(mwd(kptmw+3))
        DO 140 j = kptmw + 4, kl
          xmwa = xmwa*xbase + dble(mwa(j))
          xmwd = xmwd*xbase + dble(mwd(j))
140     CONTINUE
      ELSE
        xmwa = dble(mwa(kptmw))
        xmwd = dble(mwd(kptmw))
        DO 150 j = kptmw + 1, kl
          IF (j<=nl) THEN
            xmwa = xmwa*xbase + dble(mwa(j))
            xmwd = xmwd*xbase + dble(mwd(j))
          ELSE
            xmwa = xmwa*xbase
            xmwd = xmwd*xbase
          END IF
150     CONTINUE
      END IF
      mqdmwa = dint(xmwa*xbr)
      IF (mqdmwa<0) mqdmwa = mqdmwa - 1
      mqdmwd = dint(xmwd*xbr)
      IF (mqdmwd<0) mqdmwd = mqdmwd - 1
      mwa(kptmw) = mqdmwa
      mwa(kptmw+1) = 0
      mwa(kptmw+2) = 0
      mwd(kptmw) = mqdmwd
      mwd(kptmw+1) = 0
      mwd(kptmw+2) = 0

!             Final normalization.

      DO 160 j = kptmw - 1, 3, -1
        IF (mwa(j)<0) THEN
          mcarry = int((-mwa(j)-1)/mbase) + 1
          mwa(j) = mwa(j) + mcarry*mbase
          mwa(j-1) = mwa(j-1) - mcarry
        ELSE IF (mwa(j)>=mbase) THEN
          mcarry = -int(mwa(j)/mbase)
          mwa(j) = mwa(j) + mcarry*mbase
          mwa(j-1) = mwa(j-1) - mcarry
        END IF
        IF (mwd(j)<0) THEN
          mcarry = int((-mwd(j)-1)/mbase) + 1
          mwd(j) = mwd(j) + mcarry*mbase
          mwd(j-1) = mwd(j-1) - mcarry
        ELSE IF (mwd(j)>=mbase) THEN
          mcarry = -int(mwd(j)/mbase)
          mwd(j) = mwd(j) + mcarry*mbase
          mwd(j-1) = mwd(j-1) - mcarry
        END IF
160   CONTINUE

!             Round, affix the sign, and return.

      ma(2) = ma2
      mb(2) = mb2
      mc(1) = mc1
      mc(2) = mc2
      IF (mwa(2)==0) THEN
        mlr = 2*mwa(ndig+3) + 1
        IF (mlr>=mbase) THEN
          IF (mlr-1>mbase .AND. mwa(n1+1)<mbase-1) THEN
            IF (kround/=0 .OR. ncall>1) THEN
              mwa(n1+1) = mwa(n1+1) + 1
              mwa(n1+2) = 0
            END IF
          ELSE
            CALL fmrnd(mwa,ndig,nguard,1)
          END IF
        END IF
      ELSE
        mlr = 2*mwa(ndig+2) + 1
        IF (mlr>=mbase) THEN
          IF (mlr-1>mbase .AND. mwa(n1)<mbase-1) THEN
            IF (kround/=0 .OR. ncall>1) THEN
              mwa(n1) = mwa(n1) + 1
              mwa(n1+1) = 0
            END IF
          ELSE
            CALL fmrnd(mwa,ndig,nguard,0)
          END IF
        END IF
      END IF
      CALL fmmove(mwa,md)

      IF (mwd(2)==0) THEN
        mlr = 2*mwd(ndig+3) + 1
        IF (mlr>=mbase) THEN
          IF (mlr-1>mbase .AND. mwd(n1+1)<mbase-1) THEN
            IF (kround/=0 .OR. ncall>1) THEN
              mwd(n1+1) = mwd(n1+1) + 1
              mwd(n1+2) = 0
            END IF
          ELSE
            CALL fmrnd(mwd,ndig,nguard,1)
          END IF
        END IF
      ELSE
        mlr = 2*mwd(ndig+2) + 1
        IF (mlr>=mbase) THEN
          IF (mlr-1>mbase .AND. mwd(n1)<mbase-1) THEN
            IF (kround/=0 .OR. ncall>1) THEN
              mwd(n1) = mwd(n1) + 1
              mwd(n1+1) = 0
            END IF
          ELSE
            CALL fmrnd(mwd,ndig,nguard,0)
          END IF
        END IF
      END IF
      CALL fmmove(mwd,me)

      IF (kflag<0) THEN
        namest(ncall) = 'FMDIVD'
        CALL fmwarn
      END IF

      IF (ma2*mc2<0) md(2) = -md(2)
      IF (mb2*mc2<0) me(2) = -me(2)

      IF (kaccsw==1) THEN
        md2b = nint((ndig-1)*alogm2+log(real(abs(md(2))+1))/0.69315)
        md(0) = min(macca,maccc,md2b)
        md2b = nint((ndig-1)*alogm2+log(real(abs(me(2))+1))/0.69315)
        me(0) = min(maccb,maccc,md2b)
      ELSE
        md(0) = min(macca,maccc)
        me(0) = min(maccb,maccc)
      END IF

170   IF (ntrace/=0) THEN
        CALL fmntr(1,md,md,1)
        IF (abs(ntrace)>=1 .AND. ncall<=lvltrc) THEN
          IF (ntrace<0) THEN
            CALL fmntrj(me,ndig)
          ELSE
            CALL fmprnt(me)
          END IF
        END IF
      END IF
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmdivd
    SUBROUTINE fmdivi(ma,ival,mb)

!  MB = MA / IVAL

!  Divide FM number MA by one word integer IVAL.

!  This routine is faster than FMDIV when the divisor is less than
!  MXBASE (the square root of the largest integer).
!  When IVAL is not less than MXBASE, FMDIV2 is used.  In this case,
!  if IVAL is known to be a product of two integers less than
!  MXBASE, it is usually faster to make two calls to FMDIVI with
!  half-word factors than one call with their product.

      IMPLICIT NONE

!             Scratch array usage during FMDIVI:   M01

! .. Intrinsic Functions ..
      INTRINSIC abs, log, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: macca, md2b
! ..
! .. External Subroutines ..
      EXTERNAL fmdivn, fmntr, fmntri
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      kflag = 0
      macca = ma(0)
      ncall = ncall + 1
      IF (ntrace/=0) THEN
        namest(ncall) = 'FMDIVI'
        CALL fmntr(2,ma,ma,1)
        CALL fmntri(2,ival,0)
        CALL fmdivn(ma,ival,mb)
        IF (kaccsw==1) THEN
          md2b = nint((ndig-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
          mb(0) = min(macca,md2b)
        ELSE
          mb(0) = macca
        END IF
        CALL fmntr(1,mb,mb,1)
      ELSE
        CALL fmdivn(ma,ival,mb)
        IF (kaccsw==1) THEN
          md2b = nint((ndig-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
          mb(0) = min(macca,md2b)
        ELSE
          mb(0) = macca
        END IF
      END IF
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmdivi
    SUBROUTINE fmdivn(ma,ival,mb)

!  Internal divide by integer routine.  MB = MA / IVAL

      IMPLICIT NONE

!             Scratch array usage during FMDIVN:   M01

! .. Intrinsic Functions ..
      INTRINSIC abs, int, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma1, ma2, mkt, mlr, modint, mvalp
      INTEGER :: j, ka, kb, kl, kpt, kptwa, n1, nguard, nmval, nv2
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmdiv2, fmeq, fmim, fmmove, fmrnd, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
!             Check for special cases.

      IF (mblogs/=mbase) CALL fmcons
      n1 = ndig + 1
      IF (ma(1)==munkno .OR. ival==0) THEN
        ma1 = ma(1)
        CALL fmim(0,mb)
        mb(1) = munkno
        mb(2) = 1
        mb(0) = nint(ndig*alogm2)
        kflag = -4
        IF (ma1/=munkno) THEN
          namest(ncall) = 'FMDIVI'
          CALL fmwarn
        END IF
        RETURN
      END IF

      IF (ma(2)==0) THEN
        CALL fmeq(ma,mb)
        RETURN
      END IF

      IF (abs(ma(1))<mexpov .AND. abs(ival)>1) GO TO 20

      IF (abs(ival)==1) THEN
        DO 10 j = 0, n1
          mb(j) = ma(j)
10      CONTINUE
        mb(2) = ma(2)*ival
        IF (ma(1)==mexpov) kflag = -5
        IF (ma(1)==mexpun) kflag = -6
        RETURN
      END IF

      IF (ma(1)==mexpun) THEN
        ma2 = ma(2)
        CALL fmim(0,mb)
        mb(1) = mexpun
        mb(2) = 1
        IF ((ma2<0 .AND. ival>0) .OR. (ma2>0 .AND. ival<0)) mb(2) = -1
        mb(0) = nint(ndig*alogm2)
        kflag = -6
        RETURN
      END IF

      IF (ma(1)==mexpov) THEN
        CALL fmim(0,mb)
        mb(1) = munkno
        mb(2) = 1
        mb(0) = nint(ndig*alogm2)
        namest(ncall) = 'FMDIVI'
        kflag = -4
        CALL fmwarn
        RETURN
      END IF

!             NGUARD is the number of guard digits used.

20    IF (ncall>1) THEN
        nguard = ngrd21
      ELSE
        nguard = ngrd52
      END IF

!             If ABS(IVAL).GE.MXBASE use FMDIV.

      mvalp = abs(ival)
      nmval = int(mvalp)
      nv2 = nmval - 1
      IF (abs(ival)>mxbase .OR. nmval/=abs(ival) .OR. nv2/=abs(ival)-1) THEN
        CALL fmim(ival,m01)
        CALL fmdiv2(ma,m01,mb)
        RETURN
      END IF

!             Work with positive numbers.

      ma2 = ma(2)
      ma(2) = abs(ma(2))

!             Find the first significant digit of the quotient.

      mkt = ma(2)
      IF (mkt>=mvalp) THEN
        kpt = 2
        GO TO 50
      END IF
      DO 30 j = 3, n1
        mkt = mkt*mbase + ma(j)
        IF (mkt>=mvalp) THEN
          kpt = j
          GO TO 50
        END IF
30    CONTINUE
      kpt = n1

40    kpt = kpt + 1
      mkt = mkt*mbase
      IF (mkt<mvalp) GO TO 40

!             Do the rest of the division.

50    ka = kpt + 1
      mwa(1) = ma(1) + 2 - kpt
      mwa(2) = int(mkt/mvalp)
      modint = mkt - mwa(2)*mvalp
      kptwa = 2
      IF (ka<=n1) THEN
        kl = 3 - ka

!             (Inner Loop)

        DO 60 j = ka, n1
          mkt = modint*mbase + ma(j)
          mwa(kl+j) = int(mkt/mvalp)
          modint = mkt - mwa(kl+j)*mvalp
60      CONTINUE
        kptwa = kl + n1
      END IF

      ka = kptwa + 1
      kb = n1 + nguard
      DO 70 j = ka, kb
        mkt = modint*mbase
        mwa(j) = int(mkt/mvalp)
        modint = mkt - mwa(j)*mvalp
70    CONTINUE

!             Round the result, put the sign on MB and return.

      ma(2) = ma2
      mlr = 2*mwa(ndig+2) + 1
      IF (mlr>=mbase) THEN
        IF (mlr-1>mbase .AND. mwa(n1)<mbase-1) THEN
          IF (kround/=0 .OR. ncall>1) THEN
            mwa(n1) = mwa(n1) + 1
            mwa(n1+1) = 0
          END IF
        ELSE
          CALL fmrnd(mwa,ndig,nguard,0)
        END IF
      END IF
      CALL fmmove(mwa,mb)

      IF (kflag<0) THEN
        namest(ncall) = 'FMDIVI'
        CALL fmwarn
      END IF

      IF ((ma2<0 .AND. ival>0) .OR. (ma2>0 .AND. ival<0)) mb(2) = -mb(2)

      RETURN
    END SUBROUTINE fmdivn
    SUBROUTINE fmdm(x,ma)

!  Internal routine for converting double precision to multiple
!  precision.  Called by FMDPM.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dint, int, min, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL (KIND(0.0D0)) :: x
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: mk, mn
      REAL (KIND(0.0D0)) :: one, xbase, y, yt
      INTEGER :: j, k, n1, ne
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmdbl, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      kflag = 0
      n1 = ndig + 1

      one = 1.0D0
      xbase = mbase
      k = 0

!             NE-1 is the number of words at the current precision and
!             base roughly equal to machine precision.

      ne = int(dlogeb) + 3
      y = x
      IF (x<0.0) y = -x

      IF (x==0.0) THEN
        DO 10 j = 1, n1
          ma(j) = 0
10      CONTINUE
        GO TO 140
      END IF

!             Get the exponent.

      IF (y>one) THEN
        IF (y/xbase<y) THEN
20        k = k + 1
          y = y/xbase
          IF (y>one) GO TO 20
          IF (y<one) THEN
            ma(1) = k
            GO TO 80
          END IF
          GO TO 60
        ELSE
          DO 30 j = 1, ndig + 1
            ma(j) = 0
30        CONTINUE
          ma(1) = munkno
          ma(2) = 1
          ma(0) = nint(ndig*alogm2)
          kflag = -4
          CALL fmwarn
          RETURN
        END IF
      END IF

      IF (y<one) THEN
        IF (y*xbase>y) THEN
40        k = k - 1
          y = y*xbase
          IF (y<one) GO TO 40
          IF (y>one) THEN
            k = k + 1
            y = y/xbase
            ma(1) = k
            GO TO 80
          END IF
        ELSE
          DO 50 j = 1, ndig + 1
            ma(j) = 0
50        CONTINUE
          ma(1) = munkno
          ma(2) = 1
          ma(0) = nint(ndig*alogm2)
          kflag = -4
          CALL fmwarn
          RETURN
        END IF
      END IF

60    ma(1) = k + 1
      ma(2) = 1
      DO 70 j = 3, n1
        ma(j) = 0
70    CONTINUE
      GO TO 140

!             Build the rest of the number.

80    DO 90 j = 2, ne
        y = y*xbase
        mk = dint(y)
        yt = -mk
        CALL fmdbl(y,yt,y)
        ma(j) = mk
        IF (j>=n1) GO TO 110
90    CONTINUE
      k = ne + 1
      DO 100 j = k, n1
        ma(j) = 0
100   CONTINUE

!             Normalize.

110   IF (abs(ma(2))>=mbase) THEN
        k = n1 + 1
        DO 120 j = 3, n1
          k = k - 1
          ma(k) = ma(k-1)
120     CONTINUE
        mn = dint(ma(2)/mbase)
        ma(3) = ma(2) - mn*mbase
        ma(2) = mn
        ma(1) = ma(1) + 1
        GO TO 140
      END IF

      IF (ma(2)==0) THEN
        DO 130 j = 2, ndig
          ma(j) = ma(j+1)
130     CONTINUE
        ma(1) = ma(1) - 1
        ma(n1) = 0
      END IF

140   IF (x<0.0) ma(2) = -ma(2)
      ma(0) = min(nint((ne-1)*alogm2),nint(ndig*alogm2))
      RETURN
    END SUBROUTINE fmdm
    SUBROUTINE fmdm2(x,ma)

!  Internal routine for converting double precision to multiple
!  precision.  Called by FMDP2M.

      IMPLICIT NONE

!             Scratch array usage during FMDM2:   M01 - M04

! .. Intrinsic Functions ..
      INTRINSIC abs, dble, int, log, max, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL (KIND(0.0D0)) :: x
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: two20, y
      INTEGER :: j, jexp, k, kexp, kreslt, n1, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmi2m, fmipwr, fmmpy, fmntr, fmrslt, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
!             Increase the working precision.

      ndsave = ndig
      IF (ncall==1) THEN
        k = max(ngrd21,1)
        ndig = max(ndig+k,2)
        IF (ndig>ndg2mx) THEN
          kflag = -9
          CALL fmwarn
          ndig = ndsave
          kreslt = 12
          CALL fmrslt(ma,ma,ma,kreslt)
          IF (ntrace/=0) CALL fmntr(1,ma,ma,1)
          ncall = ncall - 1
          RETURN
        END IF
      END IF

      IF (mblogs/=mbase) CALL fmcons
      kflag = 0
      n1 = ndig + 1

      IF (x==0.0D0) THEN
        DO 10 j = 1, n1
          ma(j) = 0
10      CONTINUE
        GO TO 60
      END IF

      y = abs(x)
      two20 = 1048576.0D0

!             If this power of two is not representable at the current
!             base and precision, use a smaller one.

      IF (int(ndig*alogm2)<20) THEN
        k = int(ndig*alogm2)
        two20 = 1.0D0
        DO 20 j = 1, k
          two20 = two20*2.0D0
20      CONTINUE
      END IF

      kexp = 0
      IF (y>two20) THEN
30      y = y/two20
        kexp = kexp + 1
        IF (y>two20) GO TO 30
      ELSE IF (y<1.0D0) THEN
40      y = y*two20
        kexp = kexp - 1
        IF (y<1.0D0) GO TO 40
      END IF

      k = int(two20)
      CALL fmi2m(k,m04)
      k = int(y)
      CALL fmi2m(k,m02)
      y = (y-dble(k))*two20
      jexp = 0

50    k = int(y)
      CALL fmi2m(k,m03)
      CALL fmmpy(m02,m04,m02)
      jexp = jexp + 1
      CALL fmadd(m02,m03,m02)
      y = (y-dble(k))*two20
      IF (jexp<=1000 .AND. y/=0.0D0) GO TO 50

      k = kexp - jexp
      CALL fmipwr(m04,k,m03)
      CALL fmmpy(m02,m03,ma)

60    IF (x<0.0) ma(2) = -ma(2)
      ma(0) = nint((ndsave-1)*alogm2+log(real(abs(ma(2))+1))/0.69315)
      ndig = ndsave
      RETURN
    END SUBROUTINE fmdm2
    SUBROUTINE fmdp2m(x,ma)

!  MA = X

!  Convert a double precision floating point number to FM format.

!  This version tries to convert the double precision machine
!  number to FM with accuracy of nearly full FM precision.
!  If conversion to FM with approximately double precision accuracy
!  is good enough, FMDPM is faster and uses less scratch space.

!  This routine assumes the machine's base for double precision is
!  a power of two.

      IMPLICIT NONE

!             Scratch array usage during FMDP2M:  M01 - M04

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL (KIND(0.0D0)) :: x
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. External Subroutines ..
      EXTERNAL fmdm2, fmntr, fmntrr
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      namest(ncall) = 'FMDP2M'
      IF (ntrace/=0) CALL fmntrr(2,x,1)

      CALL fmdm2(x,ma)

      IF (ntrace/=0) CALL fmntr(1,ma,ma,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmdp2m
    SUBROUTINE fmdpm(x,ma)

!  MA = X

!  Convert a double precision floating point number to FM format.

!  In general, the relative accuracy of the FM number returned is only
!  the relative accuracy of a machine precision number.  This may be
!  true even if X can be represented exactly in the machine floating
!  point number system.

!  This version is faster than FMDP2M, but often less accurate.
!  No scratch arrays are used.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, int
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL (KIND(0.0D0)) :: x
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: y, yt
      INTEGER :: k
! ..
! .. External Subroutines ..
      EXTERNAL fmdivi, fmdm, fmim, fmntr, fmntrr
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      namest(ncall) = 'FMDPM '
      IF (ntrace/=0) CALL fmntrr(2,x,1)

!             Check to see if X is exactly a small integer.  If so,
!             converting as an integer is better.
!             Also see if X is exactly a small integer divided by
!             a small power of two.

      y = 1048576.0D0
      IF (abs(x)<y) THEN
        k = int(x)
        y = k
        IF (y==x) THEN
          CALL fmim(k,ma)
          GO TO 10
        END IF
      END IF
      IF (abs(x)<1.0D0) THEN
        y = 4096.0D0*x
        k = int(y)
        yt = k
        IF (y==yt) THEN
          CALL fmim(k,ma)
          CALL fmdivi(ma,4096,ma)
          GO TO 10
        END IF
      END IF

      CALL fmdm(x,ma)

10    IF (ntrace/=0) CALL fmntr(1,ma,ma,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmdpm
    SUBROUTINE fmentr(nroutn,ma,mb,nargs,mc,kreslt,ndsave,mxsave,kasave,kovun)

!  Do the argument checking and increasing of precision and overflow
!  threshold upon entry to an FM routine.

!  NROUTN - routine name of calling routine
!  MA     - first input argument
!  MB     - second input argument (optional)
!  NARGS  - number of input arguments
!  MC     - result argument
!  KRESLT - returned nonzero if the input arguments give the result
!           immediately (e.g., MA*0 or OVERFLOW*MB)
!  NDSAVE - saves the value of NDIG after NDIG is increased
!  MXSAVE - saves the value of MXEXP
!  KASAVE - saves the value of KACCSW
!  KOVUN  - returned nonzero if an input argument is (+ or -) overflow
!           or underflow.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC max, min
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL (KIND(0.0D0)) :: mxsave
      INTEGER :: kasave, kovun, kreslt, nargs, ndsave
      CHARACTER (6) :: nroutn
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccab
      INTEGER :: k
! ..
! .. External Subroutines ..
      EXTERNAL fmargs, fmcons, fmdivi, fmeq2, fmi2m, fmntr, fmpi, fmrslt, &
        fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      namest(ncall) = nroutn
      IF (ntrace/=0) CALL fmntr(2,ma,mb,nargs)
      CALL fmargs(nroutn,nargs,ma,mb,kreslt)

      IF (mblogs/=mbase) CALL fmcons
      kovun = 0
      IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
      IF (nargs==2) THEN
        IF (mb(1)==mexpov .OR. mb(1)==mexpun) kovun = 1
      END IF

!             Increase the working precision.

      ndsave = ndig
      IF (ncall==1) THEN
        k = max(ngrd52-1,2)
        ndig = max(ndig+k,2)
        IF (ndig>ndg2mx) THEN
          kflag = -9
          CALL fmwarn
          kreslt = 12
          ndig = ndsave
        END IF
      END IF

      IF (kreslt/=0) THEN
        maccab = ma(0)
        IF (nargs==2) maccab = min(maccab,mb(0))
        IF (kreslt==9 .OR. kreslt==10 .OR. kreslt>=13) THEN
          IF (krad==1) THEN
            CALL fmpi(mc)
          ELSE
            CALL fmi2m(180,mc)
          END IF
          IF (kreslt<=10) CALL fmdivi(mc,2,mc)
          IF (kreslt>=14) CALL fmdivi(mc,4,mc)
          CALL fmeq2(mc,mc,ndig,ndsave,1)
          ndig = ndsave
          IF (kreslt==9 .OR. kreslt==14) mc(2) = -mc(2)
          mc(0) = maccab
          IF (ntrace/=0) CALL fmntr(1,mc,mc,1)
          kasave = kaccsw
          mxsave = mxexp
          ncall = ncall - 1
          RETURN
        END IF

        ndig = ndsave
        CALL fmrslt(ma,mb,mc,kreslt)
        IF (ntrace/=0) CALL fmntr(1,mc,mc,1)
        kasave = kaccsw
        mxsave = mxexp
        ncall = ncall - 1
        RETURN
      END IF

      kasave = kaccsw
      kaccsw = 0

!             Extend the overflow/underflow threshold.

      mxsave = mxexp
      mxexp = mxexp2
      RETURN
    END SUBROUTINE fmentr
    SUBROUTINE fmeq(ma,mb)

!  MB = MA

!  This is the standard form of equality, where MA and MB both
!  have precision NDIG.  Use FMEQU for assignments that also
!  change precision.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: j
! ..
! .. External Subroutines ..
      EXTERNAL fmtrap
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      DO 10 j = 0, ndig + 1
        mb(j) = ma(j)
10    CONTINUE

!             Check for overflow or underflow.

      IF (abs(mb(1))>mxexp) THEN
        IF (mb(1)/=munkno .OR. mb(2)/=1) THEN
          ncall = ncall + 1
          CALL fmtrap(mb)
          ncall = ncall - 1
        END IF
        IF (mb(1)==munkno) kflag = -4
      END IF

      RETURN
    END SUBROUTINE fmeq
    SUBROUTINE fmeq2(ma,mb,nda,ndb,ksame)

!  Set MB (having NDB digits) equal to MA (having NDA digits).

!  If MA and MB are the same array, setting KSAME = 1 before calling
!  FMEQ2 gives faster performance.

!  If MB has less precision than MA the result is rounded to NDB digits.

!  If MB has more precision the result has zero digits padded on the
!  right.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dint, int, log, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ksame, nda, ndb
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: m2, macca, mb2, mkt
      INTEGER :: j, jt, k, kb, l, n1, ndg
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmtrap, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      macca = ma(0)

!             Check for precision in range.

      IF (nda<1 .OR. nda>ndg2mx .OR. ndb<1 .OR. ndb>ndg2mx) THEN
        ncall = ncall + 1
        namest(ncall) = 'FMEQU '
        kflag = -1
        CALL fmwarn
        WRITE (kw,90000) nda, ndb
        DO 10 j = 1, ndig + 1
          mb(j) = 0
10      CONTINUE
        mb(1) = munkno
        mb(2) = 1
        mb(0) = nint(ndig*alogm2)
        ncall = ncall - 1
        RETURN
      END IF

!             Check for special symbols.

      kflag = 0
      IF (abs(ma(1))>=mexpov) THEN
        DO 20 j = 2, ndb
          mb(j+1) = 0
20      CONTINUE
        mb(1) = ma(1)
        mb(2) = ma(2)
        GO TO 150
      END IF

      IF (ndb==nda) GO TO 100

      IF (ndb>nda) GO TO 120

!             Round to NDB digits.

      ndg = ndb
      n1 = ndb + 1
      IF (ksame/=1) THEN
        DO 30 j = 1, n1
          mb(j) = ma(j)
30      CONTINUE
      END IF
      IF (ndg<1 .OR. (kround==0 .AND. ncall<=1)) GO TO 150

      l = ndb + 2
      IF (2*(ma(l)+1)<mbase) GO TO 150
      m2 = 2
      IF (int(mbase-dint(mbase/m2)*m2)==0) THEN
        IF (2*ma(l)<mbase) GO TO 150
        IF (2*ma(l)==mbase) THEN
          IF (l<=nda) THEN
            DO 40 j = l, nda
              IF (ma(j+1)>0) GO TO 60
40          CONTINUE
          END IF

!                       Round to even.

          IF (int(mb(n1)-dint(mb(n1)/m2)*m2)==0) GO TO 150
        END IF
      ELSE
        IF (2*ma(l)+1==mbase) THEN
          IF (l<=nda) THEN
            DO 50 j = l, nda
              IF (2*(ma(j+1)+1)<mbase) GO TO 150
              IF (2*ma(j+1)>mbase) GO TO 60
50          CONTINUE
            GO TO 150
          END IF
        END IF
      END IF

60    mb(ndg+1) = mb(ndg+1) + 1
      mb(ndg+2) = 0

!             Check whether there was a carry in the rounded digit.

      mb2 = mb(2)
      mb(2) = abs(mb(2))
      kb = ndg + 1
      IF (kb>=3) THEN
        k = kb + 1
        DO 70 j = 3, kb
          k = k - 1
          IF (mb(k)<mbase) GO TO 90
          mkt = dint(mb(k)/mbase)
          mb(k-1) = mb(k-1) + mkt
          mb(k) = mb(k) - mkt*mbase
70      CONTINUE
      END IF

!             If there is a carry in the first digit then the exponent
!             must be adjusted and the number shifted right.

      IF (mb(2)<mbase) GO TO 90
      IF (kb>=4) THEN
        k = kb + 1
        DO 80 j = 4, kb
          k = k - 1
          mb(k) = mb(k-1)
80      CONTINUE
      END IF

      mkt = dint(mb(2)/mbase)
      IF (kb>=3) mb(3) = mb(2) - mkt*mbase
      mb(2) = mkt
      mb(1) = mb(1) + 1

90    IF (mb2<0) mb(2) = -mb(2)
      GO TO 150

!             MA and MB have the same precision.

100   IF (ksame/=1) THEN
        DO 110 j = 1, nda + 1
          mb(j) = ma(j)
110     CONTINUE
      END IF
      GO TO 150

!             Extend to NDB digits by padding with zeros.

120   IF (ksame/=1) THEN
        DO 130 j = 1, nda + 1
          mb(j) = ma(j)
130     CONTINUE
      END IF
      DO 140 j = nda + 2, ndb + 1
        mb(j) = 0
140   CONTINUE

!             Check for overflow or underflow.

150   IF (abs(mb(1))>mxexp) THEN
        IF (mb(1)/=munkno .OR. mb(2)/=1) THEN
          ncall = ncall + 1
          CALL fmtrap(mb)
          ncall = ncall - 1
        END IF
        IF (mb(1)==munkno) kflag = -4
      END IF

      IF (kaccsw==1) THEN
        jt = nint(log(real(abs(mb(2))+1))/0.69315)
        IF (ndb>nda) THEN
          mb(0) = nint((ndb-1)*alogm2+jt)
        ELSE
          mb(0) = min(nint((ndb-1)*alogm2+jt),int(macca))
        END IF
      ELSE
        mb(0) = ma(0)
      END IF
      RETURN
90000 FORMAT (/' The two precisions in FMEQU were NDA =',I10,' NDB =',I10/)
    END SUBROUTINE fmeq2
    SUBROUTINE fmequ(ma,mb,nda,ndb)

!  Set MB (having NDB digits) equal to MA (having NDA digits).

!  If MB has less precision than MA, the result is rounded to
!  NDB digits.

!  If MB has more precision, the result has its precision extended
!  by padding with zero digits on the right.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: nda, ndb
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. External Subroutines ..
      EXTERNAL fmeq2
! ..
      CALL fmeq2(ma,mb,nda,ndb,0)

      RETURN
    END SUBROUTINE fmequ
    SUBROUTINE fmexit(mt,mc,ndsave,mxsave,kasave,kovun)

!  Upon exit from an FM routine the result MT (having precision NDIG)
!  is rounded and returned in MC (having precision NDSAVE).
!  The values of NDIG, MXEXP, and KACCSW are restored.
!  KOVUN is nonzero if one of the routine's input arguments was overflow
!  or underflow.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL (KIND(0.0D0)) :: mxsave
      INTEGER :: kasave, kovun, ndsave
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: mc(0:lunpck), mt(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: kfsave, kwrnsv
! ..
! .. External Subroutines ..
      EXTERNAL fmeq2, fmntr, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      kwrnsv = kwarn
      kwarn = 0
      mxexp = mxsave
      kfsave = kflag
      CALL fmeq2(mt,mc,ndig,ndsave,0)
      IF (kflag/=-5 .AND. kflag/=-6) kflag = kfsave
      ndig = ndsave
      kwarn = kwrnsv
      IF (kflag==1) kflag = 0
      IF ((mc(1)==munkno .AND. kflag/=-9) .OR. (mc(1)==mexpun .AND. kovun==0) &
        .OR. (mc(1)==mexpov .AND. kovun==0)) CALL fmwarn
      IF (ntrace/=0) CALL fmntr(1,mc,mc,1)
      ncall = ncall - 1
      kaccsw = kasave
      RETURN
    END SUBROUTINE fmexit
    SUBROUTINE fmexp(ma,mb)

!  MB = EXP(MA)

      IMPLICIT NONE

!             Scratch array usage during FMEXP:   M01 - M03

! .. Intrinsic Functions ..
      INTRINSIC abs, int, log, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: m1, ma1, ma2, macca, macmax, mxsave
      REAL :: xma, xov
      INTEGER :: iextra, j, k, kasave, kovun, kreslt, kt, kwrnsv, ndmb, &
        ndsave, ndsv, nmethd
      CHARACTER (155) :: string
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmcsh2, fmdiv, fmentr, fmeq2, fmexit, fmexp2, &
        fmi2m, fmim, fmint, fmipwr, fmm2i, fmmpy, fmntr, fmrslt, fmsnh2, &
        fmsqr, fmsqrt, fmst2m, fmsub, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (abs(ma(1))>mexpab .OR. ma(2)==0) THEN
        CALL fmentr('FMEXP ',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        namest(ncall) = 'FMEXP '
        IF (ntrace/=0) CALL fmntr(2,ma,ma,1)
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
        ndsave = ndig
        IF (ncall==1) THEN
          k = max(ngrd52-1,2)
          ndig = max(ndig+k,2)
          IF (ndig>ndg2mx) THEN
            kflag = -9
            CALL fmwarn
            ndig = ndsave
            kreslt = 12
            CALL fmrslt(ma,ma,mb,kreslt)
            IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
            ncall = ncall - 1
            RETURN
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 0
        mxsave = mxexp
        mxexp = mxexp2
      END IF

      ma1 = ma(1)
      ma2 = ma(2)

      macca = ma(0)
      CALL fmeq2(ma,mb,ndsave,ndig,0)
      mb(0) = nint(ndig*alogm2)

!             Check for obvious underflow or overflow.
!             XOV is LN(LN(slightly above overflow))
!             XMA is LN(LN(EXP(MA))) approximately.

      xov = log(1.01*real(mxexp)) + log(alogmb)
      m1 = 1
      xma = log(real(max(abs(ma2),m1))) - alogmb + real(ma1)*alogmb

10    IF (xma>=xov) THEN
        CALL fmim(0,mb)
        IF (ma2>0) THEN
          kflag = -5
          mb(1) = mexpov
          mb(2) = 1
          mb(0) = nint(ndig*alogm2)
        ELSE
          kflag = -6
          mb(1) = mexpun
          mb(2) = 1
          mb(0) = nint(ndig*alogm2)
        END IF
        ndig = ndsave
        mxexp = mxsave
        kaccsw = kasave
        CALL fmwarn
        IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
        ncall = ncall - 1
        RETURN
      END IF

!             Split MA into integer and fraction parts.
!             Work with a positive argument.
!             M02 = integer part of ABS(MA)
!             MB  = fraction part of ABS(MA)

      mb(2) = abs(mb(2))
      CALL fmint(mb,m02)
      CALL fmsub(mb,m02,mb)

!             If the integer part is not zero, use FMIPWR to compute
!             E**(M02).  If M02 is too large to represent as a one word
!             integer, the definition of MXEXP insures that E**(M02)
!             overflows or underflows.

      kwrnsv = kwarn
      kwarn = 0
      CALL fmm2i(m02,kt)
      kwarn = kwrnsv
      IF (kflag/=0) THEN
        xma = xov
        GO TO 10
      END IF
      IF (kt>0) THEN

!             Compute IEXTRA, the number of extra digits required
!             to get EXP(KT) correct to the current precision.

        iextra = int(log(real(kt))/alogmb+0.5)
        IF (iextra>0 .AND. ndig+iextra<=ndg2mx) THEN
          CALL fmeq2(mb,mb,ndig,ndig+iextra,1)
        END IF
        ndig = ndig + iextra
        IF (ndig>ndg2mx) THEN
          kflag = -9
          CALL fmwarn
          mb(1) = munkno
          mb(2) = 1
          mb(0) = nint(ndig*alogm2)
          DO 20 j = 2, ndsave
            mb(j+1) = 0
20        CONTINUE
          ndig = ndig - iextra
          CALL fmexit(mb,mb,ndsave,mxsave,kasave,kovun)
          RETURN
        END IF

!             Check whether the current precision of e is large
!             enough.

        IF (mbse/=mbase .OR. ndig>ndige) THEN
          ndmb = int(150.0*2.302585/alogmb)
          IF (ndmb>=ndig) THEN
            ndsv = ndig
            ndig = min(ndmb,ndg2mx)
            string = '2.718281828459045235360287471352662497757247' // &
              '09369995957496696762772407663035354759457138217852516' // &
              '6427427466391932003059921817413596629043572900334295261'
            CALL fmst2m(string,mesav)
            mesav(0) = nint(ndig*alogm2)
            mbse = mbase
            ndige = ndig
            IF (abs(mesav(1))>10) ndige = 0
            ndig = ndsv
          ELSE
            ndsv = ndig
            ndig = min(ndig+2,ndg2mx)
            CALL fmi2m(1,mesav)
            CALL fmexp2(mesav,mesav)
            mesav(0) = nint(ndig*alogm2)
            mbse = mbase
            ndige = ndig
            IF (abs(mesav(1))>10) ndige = 0
            ndig = ndsv
          END IF
        END IF

      END IF

!             Now do the fraction part of MA and combine the results.

      kwrnsv = kwarn
      kwarn = 0
      nmethd = 1
      IF (ndig>50) nmethd = 2
      IF (mb(2)/=0 .AND. kt>0 .AND. nmethd==1) THEN
        CALL fmexp2(mb,mb)
        CALL fmipwr(mesav,kt,m03)
        CALL fmmpy(mb,m03,mb)
      ELSE IF (mb(2)/=0 .AND. kt==0 .AND. nmethd==1) THEN
        CALL fmexp2(mb,mb)
      ELSE IF (mb(2)/=0 .AND. kt>0 .AND. nmethd==2) THEN
        ndsv = ndig
        ndig = min(ndig+ngrd21,ndg2mx)
        CALL fmeq2(mb,mb,ndsv,ndig,1)
        IF (mb(1)>=0) THEN
          CALL fmcsh2(mb,mb)
          CALL fmsqr(mb,m03)
          CALL fmi2m(-1,m02)
          CALL fmadd(m03,m02,m03)
          CALL fmsqrt(m03,m03)
          CALL fmadd(mb,m03,mb)
        ELSE
          CALL fmsnh2(mb,mb)
          CALL fmsqr(mb,m03)
          CALL fmi2m(1,m02)
          CALL fmadd(m03,m02,m03)
          CALL fmsqrt(m03,m03)
          CALL fmadd(mb,m03,mb)
        END IF
        ndig = ndsv
        CALL fmipwr(mesav,kt,m03)
        CALL fmmpy(mb,m03,mb)
      ELSE IF (mb(2)/=0 .AND. kt==0 .AND. nmethd==2) THEN
        ndsv = ndig
        ndig = min(ndig+ngrd21,ndg2mx)
        CALL fmeq2(mb,mb,ndsv,ndig,1)
        IF (mb(1)>=0) THEN
          CALL fmcsh2(mb,mb)
          CALL fmsqr(mb,m03)
          CALL fmi2m(-1,m02)
          CALL fmadd(m03,m02,m03)
          CALL fmsqrt(m03,m03)
          CALL fmadd(mb,m03,mb)
        ELSE
          CALL fmsnh2(mb,mb)
          CALL fmsqr(mb,m03)
          CALL fmi2m(1,m02)
          CALL fmadd(m03,m02,m03)
          CALL fmsqrt(m03,m03)
          CALL fmadd(mb,m03,mb)
        END IF
        ndig = ndsv
      ELSE IF (mb(2)==0 .AND. kt>0) THEN
        CALL fmipwr(mesav,kt,mb)
      ELSE
        CALL fmi2m(1,mb)
      END IF

!             Invert if MA was negative.

      IF (ma2<0) THEN
        CALL fmi2m(1,m02)
        CALL fmdiv(m02,mb,mb)
      END IF
      kwarn = kwrnsv

!             Round the result and return.

      macmax = nint((ndsave-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
      mb(0) = min(mb(0),macca,macmax)
      CALL fmexit(mb,mb,ndsave,mxsave,kasave,kovun)
      RETURN
    END SUBROUTINE fmexp
    SUBROUTINE fmexp2(ma,mb)

!  MB = EXP(MA)

!  Internal exponential routine (called with 0.LT.MA.LE.1).

      IMPLICIT NONE

!             Scratch array usage during FMEXP2:   M01 - M03

!             LJSUMS = 8*(LUNPCK+1) allows for up to eight concurrent
!             sums.  Increasing this value will begin to improve the
!             speed of EXP when the base is large and precision exceeds
!             about 1,500 decimal digits.

! .. Intrinsic Functions ..
      INTRINSIC int, log, max, min, nint, real, sqrt
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL :: alog2, alogt, b, t, tj, xn
      REAL (KIND(0.0D0)) :: maxval
      INTEGER :: j, j2, k, k2, kpt, ktwo, l, l2, n2, nbig, nbot, ndsav1, &
        ndsave, nterm, ntop
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmdivi, fmeq, fmeq2, fmi2m, fmipwr, fmmpy, &
        fmmpyi, fmsqr, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mjsums(0:ljsums), &
        mlbsav(0:lunpck), mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), &
        mln4(0:lunpck), mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmsums/mjsums
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      ndsave = ndig
      IF (ma(1)==1) THEN

!             Here the special case EXP(1.0) is computed.
!             Use the direct series  e = 1/0! + 1/1! + 1/2! + ...
!             Do as much of the work as possible using small integers
!             to minimize the number of FM calls.
!             Reduce NDIG while computing each term in the
!             sum as the terms get smaller.

        t = ndig
        xn = t*alogmb/log(t)
        k = int(log(xn)/alogmb)
        ndig = max(ndig+k,2)
        IF (ndig>ndg2mx) THEN
          kflag = -9
          CALL fmwarn
          mb(1) = munkno
          mb(2) = 1
          mb(0) = nint(ndig*alogm2)
          DO 10 j = 2, ndsave
            mb(j+1) = 0
10        CONTINUE
          ndig = ndsave
          RETURN
        END IF
        ndsav1 = ndig

        CALL fmi2m(2,mb)
        CALL fmi2m(1,m02)
        j = 2
        nbig = int(mxbase)

20      ntop = 1
        nbot = j
30      IF (nbot>nbig/(j+1)) GO TO 40
        j = j + 1
        ntop = j*ntop + 1
        nbot = j*nbot
        GO TO 30

40      CALL fmdivi(m02,nbot,m02)
        IF (ntop>1) THEN
          CALL fmmpyi(m02,ntop,m03)
          ndig = ndsav1
          CALL fmadd(mb,m03,mb)
          ndig = ndsav1 - int(mb(1)-m03(1))
        ELSE
          ndig = ndsav1
          CALL fmadd(mb,m02,mb)
          ndig = ndsav1 - int(mb(1)-m02(1))
        END IF
        IF (ndig<2) ndig = 2
        IF (kflag/=1) THEN
          j = j + 1
          GO TO 20
        END IF
        ndig = ndsave
        CALL fmi2m(-1,m02)
        CALL fmadd(mb,m02,m03)
        kflag = 0
        RETURN
      END IF

!             Here is the general case.  Compute EXP(MA) where
!             0 .LT. MA .LT. 1.

!             Use the direct series
!                  EXP(X) = 1 + X + X**2/2! + X**3/3! + ...

!             The argument will be halved K2 times before the series
!             is summed.  The series will be added as J2 concurrent
!             series.  The approximately optimal values of K2 and J2
!             are now computed to try to minimize the time required.
!             N2 is the approximate number of terms of the series that
!             will be needed, and L2 guard digits will be carried.

      b = real(mbase)
      k = ngrd52
      t = max(ndig-k,2)
      alog2 = real(dlogtw)
      alogt = log(t)
      tj = 0.051*alogmb*t**0.3333 + 1.85
      j2 = int(tj)
      j2 = max(1,min(j2,ljsums/ndg2mx))
      k2 = int(1.13*sqrt(t*alogmb/tj)-0.5*alogt+4.5)

      l = int(-(real(ma(1))*alogmb+log(real(ma(2))/b+ &
        real(ma(3))/(b*b)))/alog2-0.3)
      k2 = k2 - l
      IF (l<0) l = 0
      IF (k2<0) THEN
        k2 = 0
        j2 = int(.43*sqrt(t*alogmb/(alogt+real(l)*alog2))+.33)
      END IF
      IF (j2<=1) j2 = 1

      n2 = int(t*alogmb/(alogt+real(l)*alog2))
      l2 = int(log(real(n2)+2.0**k2)/alogmb)
      ndig = ndig + l2
      IF (ndig>ndg2mx) THEN
        kflag = -9
        CALL fmwarn
        mb(1) = munkno
        mb(2) = 1
        mb(0) = nint(ndig*alogm2)
        DO 50 j = 2, ndsave
          mb(j+1) = 0
50      CONTINUE
        ndig = ndsave
        RETURN
      END IF
      ndsav1 = ndig

!             Halve the argument K2 times.

      CALL fmeq2(ma,m02,ndsave,ndig,0)
      ktwo = 1
      maxval = mxbase/2
      IF (k2>0) THEN
        DO 60 j = 1, k2
          ktwo = 2*ktwo
          IF (ktwo>maxval) THEN
            CALL fmdivi(m02,ktwo,m02)
            ktwo = 1
          END IF
60      CONTINUE
        IF (ktwo>1) CALL fmdivi(m02,ktwo,m02)
      END IF

!             Sum the series X + X**2/2! + X**3/3! + ....
!             Split into J2 concurrent sums and reduce NDIG while
!             computing each term in the sum as the terms get smaller.

      CALL fmeq(m02,mb)
      nterm = 1
      DO 70 j = 1, j2
        CALL fmdivi(mb,nterm,mb)
        nterm = nterm + 1
        kpt = (j-1)*(ndig+2)
        CALL fmeq(mb,mjsums(kpt))
70    CONTINUE
      IF (m02(1)<-ndig) GO TO 100
      CALL fmipwr(m02,j2,m03)

80    CALL fmmpy(mb,m03,mb)
      DO 90 j = 1, j2
        CALL fmdivi(mb,nterm,mb)
        kpt = (j-1)*(ndsav1+2)
        ndig = ndsav1
        CALL fmadd(mjsums(kpt),mb,mjsums(kpt))
        IF (kflag/=0) GO TO 100
        ndig = ndsav1 - int(mjsums(kpt+1)-mb(1))
        IF (ndig<2) ndig = 2
        nterm = nterm + 1
90    CONTINUE
      GO TO 80

!             Put the J2 separate sums back together.

100   kflag = 0
      kpt = (j2-1)*(ndig+2)
      CALL fmeq(mjsums(kpt),m03)
      IF (j2>=2) THEN
        DO 110 j = 2, j2
          CALL fmmpy(m02,m03,m03)
          kpt = (j2-j)*(ndig+2)
          CALL fmadd(m03,mjsums(kpt),m03)
110     CONTINUE
      END IF

!             Reverse the effect of halving the argument to
!             compute EXP(MA).

      ndig = ndsav1
      IF (k2>0) THEN
        IF (ndsave<=20) THEN
          CALL fmi2m(2,m02)
          DO 120 j = 1, k2
            CALL fmadd(m03,m02,mb)
            CALL fmmpy(mb,m03,m03)
120       CONTINUE
        ELSE
          DO 130 j = 1, k2
            CALL fmsqr(m03,mb)
            CALL fmadd(m03,m03,m02)
            CALL fmadd(mb,m02,m03)
130       CONTINUE
        END IF
      END IF
      CALL fmi2m(1,m02)
      CALL fmadd(m02,m03,mb)

      CALL fmeq2(mb,mb,ndsav1,ndsave,1)
      ndig = ndsave

      RETURN
    END SUBROUTINE fmexp2
    SUBROUTINE fmform(form,ma,string)

!  Convert an FM number (MA) to a character string base 10 (STRING)
!  using character string FORM format.

!  FORM can be one of these types:  Iw,  Fw.d,  Ew.d,  1PEw.d
!       for positive integers w,d.

!  If Iw format is used and MA is not exactly an integer, then the
!  nearest integer to MA is printed.

      IMPLICIT NONE

!             Scratch array usage during FMFORM:   M01 - M02

! .. Intrinsic Functions ..
      INTRINSIC abs, index, int, len, log10, max, min, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (*) :: form, string
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: j, jf1sav, jf2sav, jpt, k1, k2, k3, kd, ksave, kwd, kwi, &
        last, lb, lengfm, lengst, lfirst, nd, nexp
      CHARACTER (20) :: formb
! ..
! .. External Subroutines ..
      EXTERNAL fmnint, fmout
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      namest(ncall) = 'FMFORM'

      ksave = kflag
      jf1sav = jform1
      jf2sav = jform2
      string = ' '
      lengfm = len(form)
      lengst = len(string)
      kwi = 75
      kwd = 40

      IF (index(form,'I')>0 .OR. index(form,'i')>0) THEN
        k1 = max(index(form,'I'),index(form,'i')) + 1
        k2 = lengfm
        WRITE (formb,90000) k2 - k1 + 1
        IF (k2>=k1) THEN
          READ (form(k1:k2),formb) kwi
        ELSE
          kwi = lengst
        END IF
        kwi = max(1,min(kwi,lengst))
        jform1 = 2
        jform2 = 0
        kwd = kwi + 21
        IF (kwd>lmbuff) GO TO 140
        CALL fmnint(ma,m02)
        IF (m02(2)/=0) THEN
          CALL fmout(m02,cmbuff,kwd)
        ELSE
          DO 10 j = 1, kwd
            cmbuff(j) = ' '
10        CONTINUE
          cmbuff(2) = '0'
        END IF
        lfirst = 1
        last = 1
        DO 20 j = 1, kwd
          IF (cmbuff(kwd+1-j)/=' ') lfirst = kwd + 1 - j
          IF (cmbuff(j)/=' ') last = j
20      CONTINUE
        jpt = 1
        IF (last-lfirst+1>kwi) GO TO 140
        IF (last<=kwi) THEN
          DO 30 j = last, lfirst, -1
            jpt = kwi - last + j
            string(jpt:jpt) = cmbuff(j)
30        CONTINUE
          DO 40 j = 1, jpt - 1
            string(j:j) = ' '
40        CONTINUE
        ELSE
          DO 50 j = lfirst, last
            jpt = kwi - last + j
            string(jpt:jpt) = cmbuff(j)
50        CONTINUE
        END IF
      ELSE IF (index(form,'F')>0 .OR. index(form,'f')>0) THEN
        k1 = max(index(form,'F'),index(form,'f')) + 1
        k2 = index(form,'.')
        k3 = lengfm
        IF (k2>k1) THEN
          WRITE (formb,90000) k2 - k1
          READ (form(k1:k2-1),formb) kwi
        ELSE
          kwi = 50
        END IF
        IF (k3>k2) THEN
          WRITE (formb,90000) k3 - k2
          READ (form(k2+1:k3),formb) kd
        ELSE
          kd = 0
        END IF
        kwi = max(1,min(kwi,lengst))
        kd = max(0,min(kd,kwi-2))
        jform1 = 2
        jform2 = kd
        nd = int(real(ndig)*log10(real(mbase))) + 1
        IF (nd<2) nd = 2
        nexp = int(2.0*log10(real(mxbase))) + 6
        lb = max(jform2+nexp,nd+nexp)
        lb = min(lb,lmbuff)
        kwd = lb
        CALL fmout(ma,cmbuff,kwd)
        lfirst = 1
        last = 1
        DO 60 j = 1, kwd
          IF (cmbuff(kwd+1-j)/=' ') lfirst = kwd + 1 - j
          IF (cmbuff(j)/=' ') last = j
60      CONTINUE
        IF (last-lfirst+1>kwi) THEN

!             Not enough room for this F format, or FMOUT converted
!             it to E format to avoid showing no significant digits.
!             See if a shortened form will fit in E format.

          nexp = int(log10((abs(real(ma(1)))+1)*log10(real(mbase))+1)+1)
          nd = kwi - nexp - 5
          IF (nd<1) THEN
            GO TO 140
          ELSE
            jform1 = 0
            jform2 = nd
            CALL fmout(ma,cmbuff,kwi)
            lfirst = 1
            last = 1
            DO 70 j = 1, kwi
              IF (cmbuff(kwi+1-j)/=' ') lfirst = kwi + 1 - j
              IF (cmbuff(j)/=' ') last = j
70          CONTINUE
          END IF
        END IF
        jpt = 1
        IF (last<=kwi) THEN
          DO 80 j = last, lfirst, -1
            jpt = kwi - last + j
            string(jpt:jpt) = cmbuff(j)
80        CONTINUE
          DO 90 j = 1, jpt - 1
            string(j:j) = ' '
90        CONTINUE
        ELSE
          DO 100 j = lfirst, last
            jpt = kwi - last + j
            string(jpt:jpt) = cmbuff(j)
100       CONTINUE
        END IF
      ELSE IF (index(form,'1PE')>0 .OR. index(form,'1pe')>0) THEN
        k1 = max(index(form,'E'),index(form,'e')) + 1
        k2 = index(form,'.')
        k3 = lengfm
        IF (k2>k1) THEN
          WRITE (formb,90000) k2 - k1
          READ (form(k1:k2-1),formb) kwi
        ELSE
          kwi = 50
        END IF
        IF (k3>k2) THEN
          WRITE (formb,90000) k3 - k2
          READ (form(k2+1:k3),formb) kd
        ELSE
          kd = 0
        END IF
        kwi = max(1,min(kwi,lengst))
        kd = max(0,min(kd,kwi-2))
        jform1 = 1
        jform2 = kd
        IF (kwi>lmbuff) GO TO 140
        CALL fmout(ma,cmbuff,kwi)
        DO 110 j = kwi, 1, -1
          IF (j>lengst) THEN
            IF (cmbuff(j)/=' ') GO TO 140
          ELSE
            string(j:j) = cmbuff(j)
          END IF
110     CONTINUE
      ELSE IF (index(form,'E')>0 .OR. index(form,'e')>0) THEN
        k1 = max(index(form,'E'),index(form,'e')) + 1
        k2 = index(form,'.')
        k3 = lengfm
        IF (k2>k1) THEN
          WRITE (formb,90000) k2 - k1
          READ (form(k1:k2-1),formb) kwi
        ELSE
          kwi = 50
        END IF
        IF (k3>k2) THEN
          WRITE (formb,90000) k3 - k2
          READ (form(k2+1:k3),formb) kd
        ELSE
          kd = 0
        END IF
        kwi = max(1,min(kwi,lengst))
        kd = max(0,min(kd,kwi-2))
        jform1 = 0
        jform2 = kd
        IF (kwi>lmbuff) GO TO 140
        CALL fmout(ma,cmbuff,kwi)
        DO 120 j = kwi, 1, -1
          IF (j>lengst) THEN
            IF (cmbuff(j)/=' ') GO TO 140
          ELSE
            string(j:j) = cmbuff(j)
          END IF
120     CONTINUE
      ELSE
        GO TO 140
      END IF

130   kflag = ksave
      jform1 = jf1sav
      jform2 = jf2sav
      ncall = ncall - 1
      RETURN

!             Error condition.

140   kflag = -8
      DO 150 j = 1, lengst
        string(j:j) = '*'
150   CONTINUE
      GO TO 130
90000 FORMAT ('(I',I5,')')
    END SUBROUTINE fmform
    SUBROUTINE fmfprt(form,ma)

!  Print an FM number (MA) on unit KW using character
!  string FORM format.

!  FORM can be one of these types:  Iw,  Fw.d,  Ew.d,  1PEw.d
!       for positive integers w,d.

!  If Iw format is used and MA is not exactly an integer, then the
!  nearest integer to MA is printed.

      IMPLICIT NONE

!             Scratch array usage during FMFPRT:   M01 - M02

! .. Intrinsic Functions ..
      INTRINSIC abs, index, int, len, log10, max, min, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (*) :: form
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: j, jf1sav, jf2sav, jpt, k, k1, k2, k3, kd, ksave, kwd, kwi, &
        last, lb, lengfm, lfirst, nd, nexp
      CHARACTER (20) :: form2, formb
! ..
! .. External Subroutines ..
      EXTERNAL fmnint, fmout
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      namest(ncall) = 'FMFPRT'

      ksave = kflag
      jf1sav = jform1
      jf2sav = jform2
      lengfm = len(form)
      kwi = 75
      kwd = 40

      IF (index(form,'I')>0 .OR. index(form,'i')>0) THEN
        k1 = max(index(form,'I'),index(form,'i')) + 1
        k2 = lengfm
        WRITE (formb,90000) k2 - k1 + 1
        IF (k2>=k1) THEN
          READ (form(k1:k2),formb) kwi
        ELSE
          kwi = 50
        END IF
        kwi = max(1,min(kwi,lmbuff-11))
        jform1 = 2
        jform2 = 0
        kwd = kwi + 21
        CALL fmnint(ma,m02)
        IF (m02(2)/=0) THEN
          CALL fmout(m02,cmbuff,kwd)
        ELSE
          DO 10 j = 1, kwd
            cmbuff(j) = ' '
10        CONTINUE
          cmbuff(2) = '0'
        END IF
        lfirst = 1
        last = 1
        DO 20 j = 1, kwd
          IF (cmbuff(kwd+1-j)/=' ') lfirst = kwd + 1 - j
          IF (cmbuff(j)/=' ') last = j
20      CONTINUE
        jpt = 1
        IF (last-lfirst+1>kwi) GO TO 130
        IF (last<=kwi) THEN
          DO 30 j = last, lfirst, -1
            jpt = kwi - last + j
            IF (jpt/=j) cmbuff(jpt) = cmbuff(j)
30        CONTINUE
          DO 40 j = 1, jpt - 1
            cmbuff(j) = ' '
40        CONTINUE
        ELSE
          DO 50 j = lfirst, last
            jpt = kwi - last + j
            IF (jpt/=j) cmbuff(jpt) = cmbuff(j)
50        CONTINUE
        END IF
      ELSE IF (index(form,'F')>0 .OR. index(form,'f')>0) THEN
        k1 = max(index(form,'F'),index(form,'f')) + 1
        k2 = index(form(1:lengfm),'.')
        k3 = lengfm
        IF (k2>k1) THEN
          WRITE (formb,90000) k2 - k1
          READ (form(k1:k2-1),formb) kwi
        ELSE
          kwi = 50
        END IF
        IF (k3>k2) THEN
          WRITE (formb,90000) k3 - k2
          READ (form(k2+1:k3),formb) kd
        ELSE
          kd = 0
        END IF
        kwi = max(1,min(kwi,lmbuff))
        kd = max(0,min(kd,kwi-2))
        jform1 = 2
        jform2 = kd
        nd = int(real(ndig)*log10(real(mbase))) + 1
        IF (nd<2) nd = 2
        nexp = int(2.0*log10(real(mxbase))) + 6
        lb = max(jform2+nexp,nd+nexp)
        lb = min(lb,lmbuff)
        kwd = lb
        CALL fmout(ma,cmbuff,kwd)
        lfirst = 1
        last = 1
        DO 60 j = 1, kwd
          IF (cmbuff(kwd+1-j)/=' ') lfirst = kwd + 1 - j
          IF (cmbuff(j)/=' ') last = j
60      CONTINUE
        IF (last-lfirst+1>kwi) THEN

!             Not enough room for this F format, or FMOUT converted
!             it to E format to avoid showing no significant digits.
!             See if a shortened form will fit in E format.

          nexp = int(log10((abs(real(ma(1)))+1)*log10(real(mbase))+1)+1)
          nd = kwi - nexp - 5
          IF (nd<1) THEN
            GO TO 130
          ELSE
            jform1 = 0
            jform2 = nd
            CALL fmout(ma,cmbuff,kwi)
            lfirst = 1
            last = 1
            DO 70 j = 1, kwi
              IF (cmbuff(kwi+1-j)/=' ') lfirst = kwi + 1 - j
              IF (cmbuff(j)/=' ') last = j
70          CONTINUE
          END IF
        END IF
        jpt = 1
        IF (last<=kwi) THEN
          DO 80 j = last, lfirst, -1
            jpt = kwi - last + j
            IF (jpt/=j) cmbuff(jpt) = cmbuff(j)
80        CONTINUE
          DO 90 j = 1, jpt - 1
            cmbuff(j) = ' '
90        CONTINUE
        ELSE
          DO 100 j = lfirst, last
            jpt = kwi - last + j
            IF (jpt/=j) cmbuff(jpt) = cmbuff(j)
100       CONTINUE
        END IF
      ELSE IF (index(form,'1PE')>0 .OR. index(form,'1pe')>0) THEN
        k1 = max(index(form,'E'),index(form,'e')) + 1
        k2 = index(form(1:lengfm),'.')
        k3 = lengfm
        IF (k2>k1) THEN
          WRITE (formb,90000) k2 - k1
          READ (form(k1:k2-1),formb) kwi
        ELSE
          kwi = 50
        END IF
        IF (k3>k2) THEN
          WRITE (formb,90000) k3 - k2
          READ (form(k2+1:k3),formb) kd
        ELSE
          kd = 0
        END IF
        kwi = max(1,min(kwi,lmbuff))
        kd = max(0,min(kd,kwi-2))
        jform1 = 1
        jform2 = kd
        CALL fmout(ma,cmbuff,kwi)
      ELSE IF (index(form,'E')>0 .OR. index(form,'e')>0) THEN
        k1 = max(index(form,'E'),index(form,'e')) + 1
        k2 = index(form(1:lengfm),'.')
        k3 = lengfm
        IF (k2>k1) THEN
          WRITE (formb,90000) k2 - k1
          READ (form(k1:k2-1),formb) kwi
        ELSE
          kwi = 50
        END IF
        IF (k3>k2) THEN
          WRITE (formb,90000) k3 - k2
          READ (form(k2+1:k3),formb) kd
        ELSE
          kd = 0
        END IF
        kwi = max(1,min(kwi,lmbuff))
        kd = max(0,min(kd,kwi-2))
        jform1 = 0
        jform2 = kd
        CALL fmout(ma,cmbuff,kwi)
      ELSE
        GO TO 130
      END IF

110   last = kwi + 1
      WRITE (form2,90010) kswide - 7
      IF (kflag/=-8) kflag = ksave
      jform1 = jf1sav
      jform2 = jf2sav
      DO 120 j = kwi, 1, -1
        IF (cmbuff(j)/=' ' .OR. j==1) THEN
          WRITE (kw,form2) (cmbuff(k),k=1,j)
          ncall = ncall - 1
          RETURN
        END IF
120   CONTINUE
      ncall = ncall - 1
      RETURN

!             Error condition.

130   kflag = -8
      DO 140 j = 1, kwi
        cmbuff(j) = '*'
140   CONTINUE
      GO TO 110
90000 FORMAT ('(I',I5,')')
90010 FORMAT (' (6X,',I3,'A1) ')
    END SUBROUTINE fmfprt
    SUBROUTINE fmgcdi(n1,n2)

!  Find the Greatest Common Divisor of N1 and N2, and return both
!  having been divided by their GCD.  Both must be positive.

! .. Intrinsic Functions ..
      INTRINSIC max, min, mod
! ..
! .. Scalar Arguments ..
      INTEGER :: n1, n2
! ..
! .. Local Scalars ..
      INTEGER :: k1, k2, k3
! ..
      k1 = max(n1,n2)
      k2 = min(n1,n2)
10    k3 = mod(k1,k2)
      IF (k3==0) THEN
        n1 = n1/k2
        n2 = n2/k2
        RETURN
      ELSE
        k1 = k2
        k2 = k3
        GO TO 10
      END IF
    END SUBROUTINE fmgcdi
    SUBROUTINE fmi2m(ival,ma)

!  MA = IVAL

!  Convert an integer to FM format.

!  The conversion is exact if IVAL is less than MBASE**NDIG,
!  otherwise the result is an approximation.

!  This routine performs the trace printing for the conversion.
!  FMIM is used to do the arithmetic.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. External Subroutines ..
      EXTERNAL fmim, fmntr, fmntri
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (ntrace/=0) THEN
        namest(ncall) = 'FMI2M '
        CALL fmntri(2,ival,1)

        CALL fmim(ival,ma)

        CALL fmntr(1,ma,ma,1)
      ELSE
        CALL fmim(ival,ma)
      END IF
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmi2m
    SUBROUTINE fmim(ival,ma)

!  MA = IVAL.  Internal integer conversion routine.

!  The conversion is exact if IVAL is less than MBASE**NDIG.
!  Otherwise FMDM is used to get an approximation.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dint, int, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: mk, ml, mval
      REAL (KIND(0.0D0)) :: x
      INTEGER :: j, jm2, kb, kb1, n1, nmval, nv2
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmdm, fmims
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      kflag = 0
      n1 = ndig + 1

      mval = abs(ival)
      nmval = int(mval)
      nv2 = nmval - 1
      IF (abs(ival)>mxbase .OR. nmval/=abs(ival) .OR. nv2/=abs(ival)-1) THEN
        CALL fmims(ival,ma)
        GO TO 50
      END IF

!              Check for small IVAL.

      IF (mval<mbase) THEN
        DO 10 j = 3, n1
          ma(j) = 0
10      CONTINUE
        ma(2) = ival
        IF (ival==0) THEN
          ma(1) = 0
        ELSE
          ma(1) = 1
        END IF
        GO TO 50
      END IF

!             Compute and store the digits, right to left.

      ma(1) = 0
      j = ndig + 1

20    mk = dint(mval/mbase)
      ml = mval - mk*mbase
      ma(1) = ma(1) + 1
      ma(j) = ml
      IF (mk>0) THEN
        mval = mk
        j = j - 1
        IF (j>=2) GO TO 20

!             Here IVAL cannot be expressed exactly.

        x = ival
        CALL fmdm(x,ma)
        RETURN
      END IF

!             Normalize MA.

      kb = n1 - j + 2
      jm2 = j - 2
      DO 30 j = 2, kb
        ma(j) = ma(j+jm2)
30    CONTINUE
      kb1 = kb + 1
      IF (kb1<=n1) THEN
        DO 40 j = kb1, n1
          ma(j) = 0
40      CONTINUE
      END IF

      IF (ival<0) ma(2) = -ma(2)

50    ma(0) = nint(ndig*alogm2)
      RETURN
    END SUBROUTINE fmim
    SUBROUTINE fmims(ival,ma)

!  MA = IVAL.  Internal integer conversion routine.

!  This routine is called when M-variable precision is less than
!  Integer precision.  This often happens when single precision
!  is chosen for M-variables.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, int, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ml
      REAL (KIND(0.0D0)) :: x
      INTEGER :: j, jm2, kb, kb1, kbase, kmk, kval, n1
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmdm
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      kflag = 0
      n1 = ndig + 1

!              Check for small IVAL.

      kval = abs(ival)
      kbase = int(mbase)
      IF (kval<kbase) THEN
        DO 10 j = 3, n1
          ma(j) = 0
10      CONTINUE
        ma(2) = ival
        IF (ival==0) THEN
          ma(1) = 0
        ELSE
          ma(1) = 1
        END IF
        GO TO 50
      END IF

!             Compute and store the digits, right to left.

      ma(1) = 0
      j = ndig + 1

20    kmk = (kval/kbase)
      ml = kval - kmk*kbase
      ma(1) = ma(1) + 1
      ma(j) = ml
      IF (kmk>0) THEN
        kval = kmk
        j = j - 1
        IF (j>=2) GO TO 20

!             Here IVAL cannot be expressed exactly.

        x = ival
        CALL fmdm(x,ma)
        RETURN
      END IF

!             Normalize MA.

      kb = n1 - j + 2
      jm2 = j - 2
      DO 30 j = 2, kb
        ma(j) = ma(j+jm2)
30    CONTINUE
      kb1 = kb + 1
      IF (kb1<=n1) THEN
        DO 40 j = kb1, n1
          ma(j) = 0
40      CONTINUE
      END IF

      IF (ival<0) ma(2) = -ma(2)

50    ma(0) = nint(ndig*alogm2)
      RETURN
    END SUBROUTINE fmims
    SUBROUTINE fminp(line,ma,la,lb)

!  Convert an array of characters to floating point multiple precision
!  format.

!  LINE is an A1 character array of length LB to be converted
!       to FM format and returned in MA.
!  LA is a pointer telling the routine where in the array to begin
!     the conversion.  This allows more than one number to be stored
!     in an array and converted in place.
!  LB is a pointer to the last character of the field for that number.

!  The input number may be in integer or any real format.

!  KESWCH = 1  causes input to FMINP with no digits before the exponent
!              letter to be treated as if there were a leading '1'.
!              This is sometimes better for interactive input:
!              'E7' converts to 10.0**7.
!         = 0  causes a leading zero to be assumed.  This gives
!              compatibility with Fortran:
!              'E7' converts to 0.0.

!  In exponential format the 'E' may also be 'D', 'Q', or 'M'.

!  So that FMINP will convert any output from FMOUT, LINE is tested
!  to see if the input is one of the special symbols +OVERFLOW,
!  -OVERFLOW, +UNDERFLOW, -UNDERFLOW, or UNKNOWN.
!  For user input the abbreviations OVFL, UNFL, UNKN may be used.

      IMPLICIT NONE

!  Simulate a finite-state automaton to scan the input line
!  and build the number.  States of the machine:

!  1.  Initial entry to the subroutine
!  2.  Sign of the number
!  3.  Scanning digits before a decimal point
!  4.  Decimal point
!  5.  Scanning digits after a decimal point
!  6.  E, D, Q, or M -- precision indicator before the exponent
!  7.  Sign of the exponent
!  8.  Scanning exponent
!  9.  Syntax error

!  Character types recognized by the machine:

!  1.  Sign (+,-)
!  2.  Numeral (0,1,...,9)
!  3.  Decimal point (.)
!  4.  Precision indicator (E,D,Q,M)
!  5.  Illegal character for number

!  All blanks are ignored.  The analysis of the number proceeds as
!  follows:  If the simulated machine is in state JSTATE and a character
!  of type JTYPE is encountered the new state of the machine is given by
!  JTRANS(JSTATE,JTYPE).

!  In this DATA statement note the array is loaded by columns.

!          State   1  2  3  4  5  6  7  8
!  Type

! .. Intrinsic Functions ..
      INTRINSIC abs, dble, ichar, int, log10, max, min, mod, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: la, lb
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
      CHARACTER (1) :: line(lb)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: m2, mndsv1
      INTEGER :: j, jstate, k, k10pwr, kasave, kdflag, kexp, kf1, kf2, kmn, &
        kof, kpower, kpt, krsave, ksign, ksignx, kstart, kstop, ktenex, &
        ktenf1, ktenf2, ktype, kuf, kuk, kval, kwrnsv, large, n2, ndsav1, &
        ndsave
! ..
! .. Local Arrays ..
      REAL (KIND(0.0D0)) :: mlv2(0:lunpck), mlv3(0:lunpck), mlv4(0:lunpck), &
        mlv5(0:lunpck)
      INTEGER :: jtrans(8,4)
      CHARACTER (1) :: kovfl(4), kunfl(4), kunkn(4), lovfl(4), lunfl(4), &
        lunkn(4)
! ..
! .. External Subroutines ..
      EXTERNAL fmadd2, fmcons, fmdiv2, fmeq, fmeq2, fmim, fminp2, fmmi, &
        fmmpy2, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
! .. Data Statements ..
      DATA jtrans/2, 9, 9, 9, 9, 7, 9, 9, 3, 3, 3, 5, 5, 8, 8, 8, 4, 4, 4, 9, &
        9, 9, 9, 9, 6, 6, 6, 6, 6, 9, 9, 9/
      DATA kovfl/'O', 'V', 'F', 'L'/, kunfl/'U', 'N', 'F', 'L'/, kunkn/'U', &
        'N', 'K', 'N'/
      DATA lovfl/'o', 'v', 'f', 'l'/, lunfl/'u', 'n', 'f', 'l'/, lunkn/'u', &
        'n', 'k', 'n'/
! ..
!             To avoid recursion, FMINP calls only internal arithmetic
!             routines (FMADD2, FMMPY2, ...), so no trace printout is
!             done during a call to FMINP.

      IF (mblogs/=mbase) CALL fmcons
      ncall = ncall + 1
      namest(ncall) = 'FMINP '

!             Raise the call stack again, since the internal
!             routines don't.

      ncall = ncall + 1
      namest(ncall) = 'FMINP '
      ndsave = ndig
      kasave = kaccsw
      kaccsw = 0
      krsave = kround
      kround = 1
      kflag = 0

!             Check for special symbols.

      kmn = 1
      kof = 1
      kuf = 1
      kuk = 1
      DO 10 j = la, lb
        kpt = ichar(line(j))
        IF (kpt>=lhash1 .AND. kpt<=lhash2) THEN
          ktype = khasht(kpt)
          IF (ktype==2) GO TO 20
        END IF
        IF (line(j)=='-') kmn = -1
        IF (line(j)==kovfl(kof) .OR. line(j)==lovfl(kof)) THEN
          kof = kof + 1
          IF (kof==5) THEN
            CALL fmim(0,ma)
            ma(1) = mexpov
            ma(2) = kmn
            ma(0) = nint(ndig*alogm2)
            GO TO 140
          END IF
        END IF
        IF (line(j)==kunfl(kuf) .OR. line(j)==lunfl(kof)) THEN
          kuf = kuf + 1
          IF (kuf==5) THEN
            CALL fmim(0,ma)
            ma(1) = mexpun
            ma(2) = kmn
            ma(0) = nint(ndig*alogm2)
            GO TO 140
          END IF
        END IF
        IF (line(j)==kunkn(kuk) .OR. line(j)==lunkn(kof)) THEN
          kuk = kuk + 1
          IF (kuk==5) THEN
            CALL fmim(0,ma)
            ma(1) = munkno
            ma(2) = 1
            ma(0) = nint(ndig*alogm2)
            GO TO 140
          END IF
        END IF
10    CONTINUE

!             Increase the working precision.

20    IF (ncall<=2) THEN
        k = ngrd52
        ndig = max(ndig+k,2)
        IF (ndig>ndg2mx) THEN
          kflag = -9
          ncall = ncall - 1
          CALL fmwarn
          ncall = ncall + 1
          ma(1) = munkno
          ma(2) = 1
          ma(0) = nint(ndig*alogm2)
          DO 30 j = 2, ndsave
            ma(j+1) = 0
30        CONTINUE
          GO TO 140
        END IF
      END IF
      ndsav1 = ndig
      m2 = 2
      mndsv1 = ndsav1
      kstart = la
      kstop = lb
      jstate = 1
      ksign = 1
      CALL fmim(0,mlv2)
      CALL fmim(0,mlv3)
      CALL fmim(0,mlv4)
      CALL fmim(0,mlv5)

!             If MBASE is a power of ten then call FMINP2 for
!             faster input conversion.

      kpower = int(log10(dble(mbase))+0.5D0)
      IF (mbase==10**kpower) THEN
        CALL fminp2(ma,line,kstart,kstop,jtrans,kpower,mlv3,mlv4,mlv5)
        GO TO 130
      END IF

      n2 = 0
      ksignx = 1
      kf1 = 0
      kf2 = 0
      kexp = 0
      ktenf1 = 1
      ktenf2 = 1
      ktenex = 1
      k10pwr = 0

!             LARGE is a threshold used in order to do as much of the
!             conversion as possible in one-word integer arithmetic.

      large = int((intmax-10)/10)

!             KDFLAG will be 1 if any digits are found before 'E'.

      kdflag = 0

!             Scan the number.

      DO 100 j = kstart, kstop
        IF (line(j)==' ') GO TO 100
        kpt = ichar(line(j))
        IF (kpt<lhash1 .OR. kpt>lhash2) THEN
          WRITE (kw,90000) line(j), kpt, lhash1, lhash2
          ktype = 5
          kval = 0
        ELSE
          ktype = khasht(kpt)
          kval = khashv(kpt)
        END IF

        IF (ktype>=5) GO TO 150

        jstate = jtrans(jstate,ktype)

        GO TO (150,40,50,100,60,70,80,90,150) jstate

!             State 2.  Sign of the number.

40      ksign = kval
        GO TO 100

!             State 3.  Digits before a decimal point.

50      kdflag = 1
        kf1 = 10*kf1 + kval
        ktenf1 = 10*ktenf1
        IF (ktenf1>large) THEN
          IF (ktenf1/=k10pwr .AND. mlv3(2)/=0) THEN
            CALL fmim(ktenf1,ma)
            k10pwr = ktenf1
          END IF
          IF (mlv3(2)==0) THEN
            CALL fmim(kf1,mlv3)
          ELSE
            ndig = int(max(m2,min(mlv3(1)+ma(1),mndsv1)))
            CALL fmmpy2(mlv3,ma,mlv3)
            ndig = ndsav1
            CALL fmim(kf1,mlv2)
            ndig = int(max(m2,min(max(mlv3(1),mlv2(1))+1,mndsv1)))
            IF (kf1/=0) CALL fmadd2(mlv3,mlv2,mlv3)
            ndig = ndsav1
          END IF
          kf1 = 0
          ktenf1 = 1
        END IF
        GO TO 100

!             State 5.  Digits after a decimal point.

60      kdflag = 1
        n2 = n2 + 1
        kf2 = 10*kf2 + kval
        ktenf2 = 10*ktenf2
        IF (ktenf2>large) THEN
          IF (ktenf2/=k10pwr .AND. mlv4(2)/=0) THEN
            CALL fmim(ktenf2,ma)
            k10pwr = ktenf2
          END IF
          IF (mlv4(2)==0) THEN
            CALL fmim(kf2,mlv4)
          ELSE
            ndig = int(max(m2,min(mlv4(1)+ma(1),mndsv1)))
            CALL fmmpy2(mlv4,ma,mlv4)
            ndig = ndsav1
            CALL fmim(kf2,mlv2)
            ndig = int(max(m2,min(max(mlv4(1),mlv2(1))+1,mndsv1)))
            IF (kf2/=0) CALL fmadd2(mlv4,mlv2,mlv4)
            ndig = ndsav1
          END IF
          kf2 = 0
          ktenf2 = 1
        END IF
        GO TO 100

!             State 6.  Precision indicator.

70      IF (kdflag==0 .AND. keswch==1) CALL fmim(1,mlv3)
        GO TO 100

!             State 7.  Sign of the exponent.

80      ksignx = kval
        GO TO 100

!             State 8.  Digits of the exponent.

90      kexp = 10*kexp + kval
        ktenex = 10*ktenex
        IF (ktenex>large) THEN
          IF (ktenex/=k10pwr .AND. mlv5(2)/=0) THEN
            CALL fmim(ktenex,ma)
            k10pwr = ktenex
          END IF
          IF (mlv5(2)==0) THEN
            CALL fmim(kexp,mlv5)
          ELSE
            ndig = int(max(m2,min(mlv5(1)+ma(1),mndsv1)))
            CALL fmmpy2(mlv5,ma,mlv5)
            ndig = ndsav1
            CALL fmim(kexp,mlv2)
            ndig = int(max(m2,min(max(mlv5(1),mlv2(1))+1,mndsv1)))
            IF (kexp/=0) CALL fmadd2(mlv5,mlv2,mlv5)
            ndig = ndsav1
          END IF
          kexp = 0
          ktenex = 1
        END IF

100   CONTINUE

!             Form the number and return.
!             MA = KSIGN*(MLV3 + MLV4/10.0**N2)*10.0**MLV5

      IF (ktenf1>1) THEN
        IF (ktenf1/=k10pwr .AND. mlv3(2)/=0) THEN
          CALL fmim(ktenf1,ma)
          k10pwr = ktenf1
        END IF
        IF (mlv3(2)==0) THEN
          CALL fmim(kf1,mlv3)
        ELSE
          ndig = int(max(m2,min(mlv3(1)+ma(1),mndsv1)))
          CALL fmmpy2(mlv3,ma,mlv3)
          ndig = ndsav1
          CALL fmim(kf1,mlv2)
          ndig = int(max(m2,min(max(mlv3(1),mlv2(1))+1,mndsv1)))
          IF (kf1/=0) CALL fmadd2(mlv3,mlv2,mlv3)
          ndig = ndsav1
        END IF
      END IF
      IF (ktenf2>1) THEN
        IF (ktenf2/=k10pwr .AND. mlv4(2)/=0) THEN
          CALL fmim(ktenf2,ma)
          k10pwr = ktenf2
        END IF
        IF (mlv4(2)==0) THEN
          CALL fmim(kf2,mlv4)
        ELSE
          ndig = int(max(m2,min(mlv4(1)+ma(1),mndsv1)))
          CALL fmmpy2(mlv4,ma,mlv4)
          ndig = ndsav1
          CALL fmim(kf2,mlv2)
          ndig = int(max(m2,min(max(mlv4(1),mlv2(1))+1,mndsv1)))
          IF (kf2/=0) CALL fmadd2(mlv4,mlv2,mlv4)
          ndig = ndsav1
        END IF
      END IF
      IF (ktenex>1) THEN
        IF (ktenex/=k10pwr .AND. mlv5(2)/=0) THEN
          CALL fmim(ktenex,ma)
          k10pwr = ktenex
        END IF
        IF (mlv5(2)==0) THEN
          CALL fmim(kexp,mlv5)
        ELSE
          ndig = int(max(m2,min(mlv5(1)+ma(1),mndsv1)))
          CALL fmmpy2(mlv5,ma,mlv5)
          ndig = ndsav1
          CALL fmim(kexp,mlv2)
          ndig = int(max(m2,min(max(mlv5(1),mlv2(1))+1,mndsv1)))
          IF (kexp/=0) CALL fmadd2(mlv5,mlv2,mlv5)
          ndig = ndsav1
        END IF
      END IF

      IF (ksignx==-1) mlv5(2) = -mlv5(2)
      IF (mlv4(2)/=0) THEN
        CALL fmim(10,mlv2)
        k = n2
        IF (mod(k,2)==0) THEN
          CALL fmim(1,ma)
        ELSE
          CALL fmeq(mlv2,ma)
        END IF

110     k = k/2
        ndig = int(max(m2,min(2*mlv2(1),mndsv1)))
        CALL fmmpy2(mlv2,mlv2,mlv2)
        IF (mod(k,2)==1) THEN
          ndig = int(max(m2,min(mlv2(1)+ma(1),mndsv1)))
          CALL fmmpy2(mlv2,ma,ma)
        END IF
        IF (k>1) GO TO 110
        ndig = ndsav1
        CALL fmdiv2(mlv4,ma,mlv4)
      END IF
      IF (mlv5(2)/=0) THEN
        CALL fmim(10,mlv2)
        kwrnsv = kwarn
        kwarn = 0
        CALL fmmi(mlv5,kexp)
        kwarn = kwrnsv
        IF (kflag/=0) GO TO 150
        k = abs(kexp)
        IF (mod(k,2)==0) THEN
          CALL fmim(1,mlv5)
        ELSE
          CALL fmeq(mlv2,mlv5)
        END IF

120     k = k/2
        ndig = int(max(m2,min(2*mlv2(1),mndsv1)))
        CALL fmmpy2(mlv2,mlv2,mlv2)
        IF (mod(k,2)==1) THEN
          ndig = int(max(m2,min(mlv2(1)+mlv5(1),mndsv1)))
          CALL fmmpy2(mlv2,mlv5,mlv5)
        END IF
        IF (k>1) GO TO 120
        ndig = ndsav1
        IF (kexp<0) THEN
          CALL fmim(1,mlv2)
          CALL fmdiv2(mlv2,mlv5,mlv5)
        END IF
      END IF
      CALL fmadd2(mlv3,mlv4,ma)
      IF (mlv5(2)/=0) CALL fmmpy2(ma,mlv5,ma)
      IF (ksign==-1) ma(2) = -ma(2)
130   CALL fmeq2(ma,ma,ndig,ndsave,1)
      IF (ma(1)==munkno) GO TO 150

140   ndig = ndsave
      kaccsw = kasave
      kround = krsave
      IF (kflag==1) kflag = 0
      ma(0) = nint(ndig*alogm2)
      ncall = ncall - 2
      RETURN

!             Error in converting the number.

150   CALL fmim(0,ma)
      ma(1) = munkno
      ma(2) = 1
      ma(0) = nint(ndig*alogm2)
      kflag = -7
      ncall = ncall - 1
      CALL fmwarn
      ncall = ncall + 1
      GO TO 140
90000 FORMAT (/' Error in input conversion.'/ &
        ' ICHAR function was out of range for the current', &
        ' dimensions.'/' ICHAR(''',A,''') gave the value ',I12, &
        ', which is outside the currently'/' dimensioned',' bounds of (',I5, &
        ':',I5,') for variables KHASHT ','and KHASHV.'/ &
        ' Re-define the two parameters ', &
        'LHASH1 and LHASH2 so the dimensions will'/' contain', &
        ' all possible output values from ICHAR.'//)
    END SUBROUTINE fminp
    SUBROUTINE fminp2(ma,line,kstart,kstop,jtrans,kpower,mlv3,mlv4,mlv5)

!  Internal routine for input conversion for a power of ten MBASE.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC ichar, int, mod, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kpower, kstart, kstop
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mlv3(0:lunpck), mlv4(0:lunpck), mlv5(0:lunpck)
      INTEGER :: jtrans(8,4)
      CHARACTER (1) :: line(kstop)
! ..
! .. Local Scalars ..
      INTEGER :: j, jstate, kdflag, kexp, kf1, kf1dig, kf2, kf2dig, kf2pt, &
        knzdig, kpt, kshift, ksign, ksignx, ktype, kval, large
! ..
! .. External Subroutines ..
      EXTERNAL fmadd2, fmdivn, fmim, fmmpy2, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      jstate = 1
      kdflag = 0
      ksign = 1
      ksignx = 1
      kf1 = 0
      knzdig = 0
      kf1dig = 0
      kf2 = 0
      kf2dig = 0
      kf2pt = 2
      kexp = 0
      large = int(intmax/10)

!             Scan the number.

      DO 70 j = kstart, kstop
        IF (line(j)==' ') GO TO 70
        kpt = ichar(line(j))
        IF (kpt<lhash1 .OR. kpt>lhash2) THEN
          WRITE (kw,90000) line(j), kpt, lhash1, lhash2
          ktype = 5
          kval = 0
        ELSE
          ktype = khasht(kpt)
          kval = khashv(kpt)
        END IF

        IF (ktype>=5) GO TO 80

        jstate = jtrans(jstate,ktype)

        GO TO (80,10,20,70,30,40,50,60,80) jstate

!             State 2.  Sign of the number.

10      ksign = kval
        GO TO 70

!             State 3.  Digits before a decimal point.

20      kdflag = 1
        kf1 = 10*kf1 + kval
        IF (kval>0 .OR. knzdig/=0) THEN
          knzdig = 1
          kf1dig = kf1dig + 1
        END IF
        IF (kf1dig==kpower) THEN
          mlv3(1) = mlv3(1) + 1
          IF (mlv3(1)<ndig) mlv3(int(mlv3(1))+1) = kf1
          kf1 = 0
          kf1dig = 0
        END IF
        GO TO 70

!             State 5.  Digits after a decimal point.

30      kdflag = 1
        IF (kf2pt>ndig+1) GO TO 70
        kf2 = 10*kf2 + kval
        kf2dig = kf2dig + 1
        IF (kf2dig==kpower) THEN
          mlv4(kf2pt) = kf2
          IF (kf2==0 .AND. kf2pt==2) THEN
            mlv4(1) = mlv4(1) - 1
          ELSE
            kf2pt = kf2pt + 1
          END IF
          kf2 = 0
          kf2dig = 0
        END IF
        GO TO 70

!             State 6.  Precision indicator.

40      IF (kdflag==0 .AND. keswch==1) CALL fmim(1,mlv3)
        GO TO 70

!             State 7.  Sign of the exponent.

50      ksignx = kval
        GO TO 70

!             State 8.  Digits of the exponent.

60      IF (kexp>=large) THEN
          IF (mlv3(2)==0 .AND. mlv4(2)==0) THEN
            CALL fmim(0,ma)
            RETURN
          END IF
          CALL fmim(0,ma)
          IF (ksignx==1) THEN
            ma(1) = mexpov
            kflag = -4
          ELSE
            ma(1) = mexpun
            kflag = -4
          END IF
          ma(2) = ksign
          ma(0) = nint(ndig*alogm2)
          ncall = ncall - 1
          CALL fmwarn
          ncall = ncall + 1
          RETURN
        END IF

        kexp = 10*kexp + kval

70    CONTINUE

!             Form the number and return.
!             MA = KSIGN*(MLV3 + MLV4)*10.0**(KSIGNX*KEXP)

      IF (kf1dig/=0) THEN
        mlv3(1) = mlv3(1) + 1
        kshift = 10**(kpower-kf1dig)
        IF (mlv3(1)<ndig) mlv3(int(mlv3(1))+1) = kf1*kshift
        IF (kshift>1) THEN
          CALL fmdivn(mlv3,kshift,mlv3)
        END IF
      END IF

      IF (kf2dig/=0) THEN
        kshift = 10**(kpower-kf2dig)
        mlv4(kf2pt) = kf2*kshift
      END IF
      IF (mlv4(2)==0) mlv4(1) = 0

      IF (kexp/=0) THEN
        IF (ksignx==1) THEN
          mlv5(1) = int(kexp/kpower) + 1
          mlv5(2) = 10**(mod(kexp,kpower))
        ELSE
          mlv5(1) = -int((kexp-1)/kpower)
          kshift = 10**(mod(kexp,kpower))
          IF (kshift>1) THEN
            mlv5(2) = mbase/kshift
          ELSE
            mlv5(2) = 1
          END IF
        END IF
      END IF

      CALL fmadd2(mlv3,mlv4,ma)
      IF (kexp>0) CALL fmmpy2(ma,mlv5,ma)
      ma(2) = ksign*ma(2)

      RETURN

!             Error in converting the number.

80    CALL fmim(0,ma)
      ma(1) = munkno
      ma(2) = 1
      ma(0) = nint(ndig*alogm2)
      RETURN
90000 FORMAT (/' Error in input conversion.'/ &
        ' ICHAR function was out of range for the current', &
        ' dimensions.'/' ICHAR(''',A,''') gave the value ',I12, &
        ', which is outside the currently'/' dimensioned',' bounds of (',I5, &
        ':',I5,') for variables KHASHT ','and KHASHV.'/ &
        ' Re-define the two parameters ', &
        'LHASH1 and LHASH2 so the dimensions will'/' contain', &
        ' all possible output values from ICHAR.'//)
    END SUBROUTINE fminp2
    SUBROUTINE fmint(ma,mb)

!  MB = INT(MA)

!  The integer part of MA is computed and returned in MB as a multiple
!  precision floating point number.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, int, log, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: macca, macmax
      INTEGER :: j, ka, kb, kreslt, n1
! ..
! .. External Subroutines ..
      EXTERNAL fmargs, fmcons, fmntr, fmrslt
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      macca = ma(0)
      kflag = 0
      ncall = ncall + 1
      namest(ncall) = 'FMINT '
      IF (ntrace/=0) CALL fmntr(2,ma,ma,1)
      IF (abs(ma(1))>mexpab) THEN
        CALL fmargs('FMINT ',1,ma,mb,kreslt)
        IF (kreslt/=0) THEN
          CALL fmrslt(ma,ma,mb,kreslt)
          IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
          ncall = ncall - 1
          RETURN
        END IF
      END IF

      n1 = ndig + 1

!             If MA is less than one in magnitude, return zero.

      IF (ma(1)<=0) THEN
        DO 10 j = 1, n1
          mb(j) = 0
10      CONTINUE
        GO TO 50
      END IF

!             If the radix point is off the right end of MA then MA is
!             already an integer.  Return MA.

      IF (ma(1)>=ndig) THEN
        DO 20 j = 1, n1
          mb(j) = ma(j)
20      CONTINUE
        GO TO 50
      END IF

!             Here MA has both integer and fraction parts.  Replace
!             the digits right of the radix point by zeros.

      ka = int(ma(1)) + 2
      kb = ka - 1
      DO 30 j = 1, kb
        mb(j) = ma(j)
30    CONTINUE

      DO 40 j = ka, n1
        mb(j) = 0
40    CONTINUE

50    IF (kaccsw==1) THEN
        macmax = nint((ndig-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
        mb(0) = min(macca,macmax)
      ELSE
        mb(0) = macca
      END IF
      IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmint
    SUBROUTINE fmipwr(ma,ival,mb)

!  MB = MA ** IVAL

!  Raise an FM number to an integer power.
!  The binary multiplication method used requires an average of
!  1.5 * LOG2(IVAL) multiplications.  MA may be negative.

      IMPLICIT NONE

!             Scratch array usage during FMIPWR:   M01

! .. Intrinsic Functions ..
      INTRINSIC abs, int, log, max, min, mod, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, macca, macmax
      REAL :: xval
      INTEGER :: j, jsign, k, kwrnsv, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmdiv, fmeq, fmeq2, fmi2m, fmim, fmmpy, fmntr, fmntri, &
        fmsqr, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      kflag = 0
      ncall = ncall + 1
      namest(ncall) = 'FMIPWR'
      IF (ntrace/=0) THEN
        CALL fmntr(2,ma,ma,1)
        CALL fmntri(2,ival,0)
      END IF

!             Check for special cases.

      IF (ma(1)==munkno .OR. (ival<=0 .AND. ma(2)==0)) THEN
        ma2 = ma(2)
        CALL fmim(0,mb)
        mb(1) = munkno
        mb(2) = 1
        mb(0) = nint(ndig*alogm2)
        kflag = -4
        IF (ival<=0 .AND. ma2==0) CALL fmwarn
        IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
        ncall = ncall - 1
        RETURN
      END IF

      IF (ival==0) THEN
        CALL fmim(1,mb)
        IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
        ncall = ncall - 1
        RETURN
      END IF

      IF (abs(ival)==1) THEN
        kwrnsv = kwarn
        kwarn = 0
        IF (ival==1) THEN
          CALL fmeq(ma,mb)
        ELSE
          CALL fmim(1,m01)
          CALL fmdiv(m01,ma,mb)
        END IF
        IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
        ncall = ncall - 1
        kwarn = kwrnsv
        RETURN
      END IF

      IF (ma(2)==0) THEN
        CALL fmeq(ma,mb)
        IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
        ncall = ncall - 1
        RETURN
      END IF

      IF (ma(1)==mexpov) THEN
        jsign = 1
        IF (ma(2)<0) jsign = -1
        CALL fmim(0,mb)
        IF (ival>0) THEN
          mb(1) = mexpov
          mb(2) = jsign**mod(ival,2)
          mb(0) = nint(ndig*alogm2)
          kflag = -5
        ELSE
          mb(1) = mexpun
          mb(2) = jsign**mod(ival,2)
          mb(0) = nint(ndig*alogm2)
          kflag = -6
        END IF
        IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
        ncall = ncall - 1
        RETURN
      END IF

      IF (ma(1)==mexpun) THEN
        jsign = 1
        IF (ma(2)<0) jsign = -1
        CALL fmim(0,mb)
        IF (ival>0) THEN
          mb(1) = mexpun
          mb(2) = jsign**mod(ival,2)
          mb(0) = nint(ndig*alogm2)
          kflag = -6
        ELSE
          mb(1) = mexpov
          mb(2) = jsign**mod(ival,2)
          mb(0) = nint(ndig*alogm2)
          kflag = -5
        END IF
        IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
        ncall = ncall - 1
        RETURN
      END IF

!             Increase the working precision.

      ndsave = ndig
      IF (ncall==1) THEN
        xval = abs(ival)
        k = int((5.0*real(dlogtn)+log(xval))/alogmb+2.0)
        ndig = max(ndig+k,2)
      ELSE
        xval = abs(ival)
        IF (xval>10.0 .OR. real(mbase)<=999.0) THEN
          k = int(log(xval)/alogmb+1.0)
          ndig = ndig + k
        END IF
      END IF
      IF (ndig>ndg2mx) THEN
        kflag = -9
        CALL fmwarn
        mb(1) = munkno
        mb(2) = 1
        mb(0) = nint(ndig*alogm2)
        DO 10 j = 2, ndsave
          mb(j+1) = 0
10      CONTINUE
        ndig = ndsave
        IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
        ncall = ncall - 1
        RETURN
      END IF

!             Initialize.

      kwrnsv = kwarn
      kwarn = 0
      k = abs(ival)

      macca = ma(0)
      CALL fmeq2(ma,m01,ndsave,ndig,0)
      m01(0) = nint(ndig*alogm2)

      IF (mod(k,2)==0) THEN
        CALL fmi2m(1,mb)
      ELSE
        CALL fmeq(m01,mb)
      END IF

!             This is the multiplication loop.

20    k = k/2
      CALL fmsqr(m01,m01)
      IF (mod(k,2)==1) CALL fmmpy(m01,mb,mb)
      IF (k>1) GO TO 20

!             Invert if the exponent is negative.

      IF (ival<0) THEN
        CALL fmi2m(1,m01)
        CALL fmdiv(m01,mb,mb)
      END IF
      kwarn = kwrnsv

!             Round the result and return.

      CALL fmeq2(mb,mb,ndig,ndsave,1)
      ndig = ndsave
      IF (kaccsw==1) THEN
        macmax = nint((ndsave-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
        mb(0) = min(mb(0),macca,macmax)
      ELSE
        mb(0) = macca
      END IF
      IF (kflag<0) CALL fmwarn
      IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmipwr
    SUBROUTINE fmlg10(ma,mb)

!  MB = LOG10(MA)

      IMPLICIT NONE

!             Scratch array usage during FMLG10:   M01 - M05

! .. Intrinsic Functions ..
      INTRINSIC abs, log, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: macca, macmax, mxsave
      INTEGER :: k, kasave, kovun, kreslt, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmdiv, fmentr, fmeq2, fmexit, fmln, fmlni, &
        fmntr, fmrslt, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (abs(ma(1))>mexpab .OR. ma(2)<=0) THEN
        CALL fmentr('FMLG10',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        namest(ncall) = 'FMLG10'
        IF (ntrace/=0) CALL fmntr(2,ma,ma,1)
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
        ndsave = ndig
        IF (ncall==1) THEN
          k = max(ngrd52-1,2)
          ndig = max(ndig+k,2)
          IF (ndig>ndg2mx) THEN
            kflag = -9
            CALL fmwarn
            ndig = ndsave
            kreslt = 12
            CALL fmrslt(ma,ma,mb,kreslt)
            IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
            ncall = ncall - 1
            RETURN
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 0
        mxsave = mxexp
        mxexp = mxexp2
      END IF

      macca = ma(0)
      CALL fmeq2(ma,mb,ndsave,ndig,0)
      mb(0) = nint(ndig*alogm2)

      CALL fmln(mb,mb)
      IF (mbase/=mbsli .OR. ndig>ndigli) THEN
        CALL fmlni(10,m03)
      ELSE
        CALL fmadd(mln1,mln3,m03)
      END IF
      CALL fmdiv(mb,m03,mb)

!             Round the result and return.

      macmax = nint((ndsave-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
      mb(0) = min(mb(0),macca,macmax)
      CALL fmexit(mb,mb,ndsave,mxsave,kasave,kovun)
      RETURN
    END SUBROUTINE fmlg10
    SUBROUTINE fmln(ma,mb)

!  MB = LOG(MA)     (Natural logarithm)

      IMPLICIT NONE

!             Scratch array usage during FMLN:   M01 - M05

! .. Intrinsic Functions ..
      INTRINSIC abs, int, log, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma1, macca, macmax, mxsave
      REAL (KIND(0.0D0)) :: y
      REAL :: x
      INTEGER :: iextra, ival, j, k, k2, k2exp, kasave, kbot, km1, kovun, &
        kreslt, kscale, kst, kwrnsv, last, n1, n3, ndsav1, ndsave, ndsv
! ..
! .. Local Arrays ..
      INTEGER :: nstack(19)
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmdig, fmdiv, fmdivi, fmdpm, fmentr, fmeq, &
        fmeq2, fmexit, fmexp, fmi2m, fmlni, fmm2dp, fmm2i, fmmpy, fmmpyi, &
        fmntr, fmrslt, fmsqr, fmsub, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (abs(ma(1))>mexpab .OR. ma(2)<=0) THEN
        CALL fmentr('FMLN  ',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        namest(ncall) = 'FMLN  '
        IF (ntrace/=0) CALL fmntr(2,ma,ma,1)
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
        ndsave = ndig
        IF (ncall==1) THEN
          k = max(ngrd52-1,2)
          ndig = max(ndig+k,2)
          IF (ndig>ndg2mx) THEN
            kflag = -9
            CALL fmwarn
            ndig = ndsave
            kreslt = 12
            CALL fmrslt(ma,ma,mb,kreslt)
            IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
            ncall = ncall - 1
            RETURN
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 0
        mxsave = mxexp
        mxexp = mxexp2
      END IF

!             If MA is close to 1, use the Taylor series:
!                   LN(1+X) = X - X**2/2 + X**3/3 - ...
!             This is faster for small X and avoids cancellation error.

!             This method is faster for moderate sized NDIG, but is
!             asymptotically slower by a factor of NDIG**(2/3) than
!             using Newton and FMEXP.  For MBASE=10,000 the Taylor
!             series is faster for NDIG less than about 150 (and is
!             used only when MA is between .9999 and 1.0001).

      IF (ma(1)==0 .OR. ma(1)==1) THEN
        x = real(mbase)
        x = x**(int(ma(1))-1)*(real(ma(2))+real(ma(3))/x)
      ELSE
        x = 2.0
      END IF
      IF (x>0.9999 .AND. x<=1.0001) THEN
        macca = ma(0)
        CALL fmeq2(ma,m03,ndsave,ndig,0)
        m03(0) = nint(ndig*alogm2)

        CALL fmi2m(-1,m01)
        CALL fmadd(m03,m01,m03)

!             The sum will be done as two concurrent series.

        ndsav1 = ndig
        CALL fmeq(m03,m04)
        CALL fmdivi(m03,2,m05)
        CALL fmsqr(m03,mb)
        CALL fmeq(m03,m02)
        kbot = 2

10      kbot = kbot + 1
        CALL fmmpy(m02,mb,m02)
        CALL fmdivi(m02,kbot,m01)
        ndig = ndsav1
        CALL fmadd(m04,m01,m04)
        ndig = max(2,ndsav1-int(m04(1)-m01(1)))
        kbot = kbot + 1
        CALL fmdivi(m02,kbot,m01)
        ndig = ndsav1
        CALL fmadd(m05,m01,m05)
        ndig = max(2,ndsav1-int(m04(1)-m01(1)))
        IF (kflag/=1) GO TO 10

        ndig = ndsav1
        CALL fmmpy(m05,m03,m05)
        CALL fmsub(m04,m05,mb)
        GO TO 70
      END IF

      ma1 = ma(1)
      macca = ma(0)
      CALL fmeq2(ma,m05,ndsave,ndig,0)
      m05(0) = nint(ndig*alogm2)

!             Compute IEXTRA, the number of extra digits required.

      CALL fmi2m(1,m04)
      CALL fmsub(m04,m05,m04)
      iextra = max(0-int(m04(1)),0)
      IF (iextra>0 .AND. ndig+iextra<=ndg2mx) THEN
        CALL fmeq2(m05,m05,ndig,ndig+iextra,1)
      END IF
      ndig = ndig + iextra
      IF (ndig>ndg2mx) THEN
        kflag = -9
        CALL fmwarn
        mb(1) = munkno
        mb(2) = 1
        mb(0) = nint(ndig*alogm2)
        DO 20 j = 2, ndsave
          mb(j+1) = 0
20      CONTINUE
        ndig = ndig - iextra
        CALL fmexit(mb,mb,ndsave,mxsave,kasave,kovun)
        RETURN
      END IF

!             Check to see if the argument is a small integer.
!             If so use FMLNI.

      km1 = 0

      kwrnsv = kwarn
      kwarn = 0
      CALL fmm2i(m05,ival)
      kwarn = kwrnsv
      IF (kflag==0 .AND. ival<mxbase) THEN
        CALL fmlni(ival,mb)
        GO TO 70
      END IF

!             See if the argument can be scaled to a small integer.

      n3 = ndig + 3
      n1 = ndig + 1
      DO 30 j = 2, n1
        IF (m05(n3-j)/=0) THEN
          last = n3 - j - 1
          GO TO 40
        END IF
30    CONTINUE

40    kscale = int(ma1) - last
      m05(1) = last
      kwrnsv = kwarn
      kwarn = 0
      CALL fmm2i(m05,ival)
      kwarn = kwrnsv
      IF (kflag==0 .AND. ival<mxbase) THEN
        CALL fmlni(ival,m04)
        IF (ival==1) km1 = 1
        k2exp = 0
        GO TO 60
      END IF

!             For the non-integer case, scale the argument to lie
!             between e/2 and e to speed up the calls to FMEXP.

      m05(1) = 1
      kscale = int(ma1) - 1
      CALL fmm2dp(m05,y)
      k2exp = int(log(2.0*real(y)/2.71828)/0.693147)
      IF (y<1.359141) THEN
        k2exp = -1
        CALL fmadd(m05,m05,m05)
        y = 2.0D0*y
      ELSE
        k2 = 2**k2exp
        CALL fmdivi(m05,k2,m05)
        y = y/k2
      END IF

!             Generate the initial approximation.

      y = log(y)
      CALL fmdpm(y,m04)
      CALL fmdig(nstack,kst)

!             Newton iteration.

      DO 50 j = 1, kst
        ndig = nstack(j)
        CALL fmexp(m04,mb)
        CALL fmsub(m05,mb,m02)
        CALL fmdiv(m02,mb,mb)
        CALL fmadd(m04,mb,m04)
50    CONTINUE
      m04(0) = nint(ndig*alogm2)

!             Compute LN(MBASE**KSCALE).

60    IF ((mbslb/=mbase .OR. ndiglb<ndig) .AND. kscale/=0) THEN
        ndsv = ndig
        ndig = min(ndig+2,ndg2mx)
        CALL fmlni(int(mbase),mlbsav)
        mbslb = mbase
        ndiglb = ndig
        IF (abs(mlbsav(1))>10) ndiglb = 0
        ndig = ndsv
      END IF

      IF (kscale/=0 .AND. km1==0) THEN
        CALL fmmpyi(mlbsav,kscale,mb)
        CALL fmadd(m04,mb,mb)
      ELSE IF (kscale/=0 .AND. km1==1) THEN
        CALL fmmpyi(mlbsav,kscale,mb)
      ELSE IF (kscale==0 .AND. km1==0) THEN
        CALL fmeq(m04,mb)
      ELSE IF (kscale==0 .AND. km1==1) THEN
        CALL fmi2m(0,mb)
      END IF

      IF (k2exp/=0) THEN
        IF (mbase/=mbsli .OR. ndig>ndigli) THEN
          CALL fmlni(2,m04)
        END IF
        CALL fmmpyi(mln1,k2exp,m04)
        CALL fmadd(mb,m04,mb)
      END IF

!             Round the result and return.

70    macmax = nint((ndsave-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
      mb(0) = min(mb(0),macca,macmax)
      CALL fmexit(mb,mb,ndsave,mxsave,kasave,kovun)
      RETURN
    END SUBROUTINE fmln
    SUBROUTINE fmlni(ival,ma)

!  MA = LOG(IVAL)

!  Compute the natural logarithm of an integer IVAL.

!  If IVAL has only powers of 2, 3, 5, and 7 in its factorization then
!  FMLNI is faster than FMLN.  Otherwise, if IVAL.GE.MXBASE (i.e., IVAL
!  does not fit in 1/2 word) then FMLN is usually faster.

!  Use FMLN instead of FMLNI if 10*IVAL would cause integer overflow.

      IMPLICIT NONE

!             Scratch array usage during FMLNI:   M01 - M03

! .. Intrinsic Functions ..
      INTRINSIC abs, dble, int, log, max, min, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      REAL :: xval
      INTEGER :: int2, j, j2, j3, j5, j7, jtemp2, jtemp3, jtemp5, jtemp7, k, &
        k2, k3, k5, k7, kasave, kdelta, last, nd, ndmb, ndsave, ndsv, nt
      CHARACTER (155) :: string
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmdivi, fmeq2, fmi2m, fmim, fmlni2, fmmpyi, &
        fmntr, fmntri, fmst2m, fmsub, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      kflag = 0
      ncall = ncall + 1
      namest(ncall) = 'FMLNI '
      IF (ntrace/=0) CALL fmntri(2,ival,1)

!             Check for special cases.

      IF (ival<=0) THEN
        CALL fmim(0,ma)
        ma(1) = munkno
        ma(2) = 1
        ma(0) = nint(ndig*alogm2)
        kflag = -4
        CALL fmwarn
        IF (ntrace/=0) CALL fmntr(1,ma,ma,1)
        ncall = ncall - 1
        RETURN
      END IF

      IF (ival==1) THEN
        CALL fmi2m(0,ma)
        IF (ntrace/=0) CALL fmntr(1,ma,ma,1)
        ncall = ncall - 1
        RETURN
      END IF

!             Increase the working precision.

      ndsave = ndig
      IF (ncall==1) THEN
        k = ngrd52
        ndig = max(ndig+k,2)
        IF (ndig>ndg2mx) THEN
          kflag = -9
          CALL fmwarn
          ma(1) = munkno
          ma(2) = 1
          ma(0) = nint(ndsave*alogm2)
          DO 10 j = 2, ndsave
            ma(j+1) = 0
10        CONTINUE
          ndig = ndsave
          IF (ntrace/=0) CALL fmntr(1,ma,ma,1)
          ncall = ncall - 1
          RETURN
        END IF
      END IF
      kasave = kaccsw
      kaccsw = 0

!             Find integers K2, K3, K5, and K7 such that
!                NT = 2**K2 * 3**K3 * 5**K5 * 7**K7
!             is a good approximation of IVAL.
!             KDELTA = ABS(IVAL - NT).

      int2 = ival
      IF (ival>intmax/100) int2 = ival/100
      kdelta = int2
      nt = 0
      k2 = 0
      k3 = 0
      k5 = 0
      k7 = 0

!             Start the search loop.

      xval = int2
      last = int(log(dble(xval))/dlogtw+2.0D0)

      jtemp7 = 1
      DO 80 j7 = 1, last
        IF (jtemp7>int2 .AND. abs(jtemp7-int2)>kdelta) GO TO 90

        jtemp5 = jtemp7
        DO 60 j5 = 1, last
          IF (jtemp5>int2 .AND. abs(jtemp5-int2)>kdelta) GO TO 70

          jtemp3 = jtemp5
          DO 40 j3 = 1, last
            IF (jtemp3>int2 .AND. abs(jtemp3-int2)>kdelta) GO TO 50

            jtemp2 = jtemp3
            DO 20 j2 = 1, last
              IF (abs(jtemp2-int2)<=kdelta) THEN
                IF (abs(jtemp2-int2)==kdelta .AND. jtemp2<int2) GO TO 30
                kdelta = abs(jtemp2-int2)
                nt = jtemp2
                k2 = j2 - 1
                k3 = j3 - 1
                k5 = j5 - 1
                k7 = j7 - 1
                IF (kdelta==0) GO TO 90
              END IF
              IF (jtemp2>int2) GO TO 30

              jtemp2 = 2*jtemp2
20          CONTINUE

30          jtemp3 = 3*jtemp3
40        CONTINUE

50        jtemp5 = 5*jtemp5
60      CONTINUE

70      jtemp7 = 7*jtemp7
80    CONTINUE

!             If IVAL was too close to the integer overflow limit,
!             restore NT to an approximation of IVAL.

90    IF (int2/=ival) THEN
        IF (nt<=int2) THEN
          nt = nt*100
          k2 = k2 + 2
          k5 = k5 + 2
        ELSE IF (nt<=ival/98) THEN
          nt = nt*98
          k2 = k2 + 1
          k7 = k7 + 2
        ELSE
          nt = nt*70
          k2 = k2 + 1
          k5 = k5 + 1
          k7 = k7 + 1
        END IF
      END IF

!             End of the search.  Now compute LN(NT) as a linear
!             combination of LN(2), LN(3), LN(5), and LN(7).

      IF (mbase/=mbsli .OR. ndig>ndigli) THEN
        ndmb = int(150.0*2.302585/alogmb)
        IF (ndmb>=ndig) THEN
          ndsv = ndig
          ndig = min(ndmb,ndg2mx)
          string = '0.693147180559945309417232121458176568075500' // &
            '13436025525412068000949339362196969471560586332699641' // &
            '8687542001481020570685733685520235758130557032670751635'
          CALL fmst2m(string,mln1)
          string = '1.098612288668109691395245236922525704647490' // &
            '55782274945173469433363749429321860896687361575481373' // &
            '2088787970029065957865742368004225930519821052801870767'
          CALL fmst2m(string,mln2)
          string = '1.609437912434100374600759333226187639525601' // &
            '35426851772191264789147417898770765776463013387809317' // &
            '9610799966303021715562899724005229324676199633616617464'
          CALL fmst2m(string,mln3)
          string = '1.945910149055313305105352743443179729637084' // &
            '72958186118845939014993757986275206926778765849858787' // &
            '1526993061694205851140911723752257677786843148958095164'
          CALL fmst2m(string,mln4)
          mbsli = mbase
          ndigli = ndig
          IF (abs(mln1(1))>10 .OR. abs(mln2(1))>10 .OR. abs(mln3( &
            1))>10 .OR. abs(mln4(1))>10) ndigli = 0
        ELSE
          ndsv = ndig
          ndig = min(ndig+2,ndg2mx)
          mbsli = mbase
          ndigli = ndig

          CALL fmlni2(1,126,mln1)
          CALL fmlni2(1,225,mln2)
          CALL fmlni2(1,2401,mln3)
          CALL fmlni2(1,4375,mln4)

!                Get Ln(2).

          CALL fmmpyi(mln1,-72,mln1)
          CALL fmmpyi(mln2,-27,ma)
          CALL fmadd(mln1,ma,mln1)
          CALL fmmpyi(mln3,19,ma)
          CALL fmadd(mln1,ma,mln1)
          CALL fmmpyi(mln4,-31,ma)
          CALL fmadd(mln1,ma,mln1)

!                Get Ln(3).

          CALL fmmpyi(mln2,-3,mln2)
          CALL fmmpyi(mln1,19,ma)
          CALL fmadd(mln2,ma,mln2)
          CALL fmsub(mln2,mln3,mln2)
          CALL fmadd(mln2,mln4,mln2)
          CALL fmdivi(mln2,12,mln2)

!                Get Ln(5).

          CALL fmsub(mln3,mln1,mln3)
          CALL fmmpyi(mln2,27,ma)
          CALL fmadd(mln3,ma,mln3)
          CALL fmmpyi(mln4,-4,ma)
          CALL fmadd(mln3,ma,mln3)
          CALL fmdivi(mln3,18,mln3)

!                Get Ln(7).

          CALL fmsub(mln1,mln4,mln4)
          CALL fmmpyi(mln2,7,ma)
          CALL fmadd(mln4,ma,mln4)
          CALL fmmpyi(mln3,-4,ma)
          CALL fmadd(mln4,ma,mln4)
        END IF
        mln1(0) = nint(ndig*alogm2)
        mln2(0) = mln1(0)
        mln3(0) = mln1(0)
        mln4(0) = mln1(0)
        IF (abs(mln1(1))>10 .OR. abs(mln2(1))>10 .OR. abs(mln3( &
          1))>10 .OR. abs(mln4(1))>10) ndigli = 0
        ndig = ndsv
      END IF

!             If NT.NE.IVAL then the final step is to compute
!             LN(IVAL/NT) and then use LN(IVAL) = LN(IVAL/NT) + LN(NT).

      IF (nt/=ival) THEN
        nd = nt - ival
        CALL fmlni2(nd,nt,ma)
      END IF

      CALL fmmpyi(mln1,k2,m02)
      CALL fmmpyi(mln2,k3,m01)
      CALL fmadd(m02,m01,m02)
      CALL fmmpyi(mln3,k5,m01)
      CALL fmadd(m02,m01,m02)
      CALL fmmpyi(mln4,k7,m01)
      IF (nt/=ival) CALL fmadd(m02,ma,m02)
      CALL fmadd(m02,m01,ma)

!             Round and move the result to MA.

      kaccsw = kasave
      CALL fmeq2(ma,ma,ndig,ndsave,1)
      ndig = ndsave
      IF (ntrace/=0) CALL fmntr(1,ma,ma,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmlni
    SUBROUTINE fmlni2(int1,int2,ma)

!  MA = LN(1 - INT1/INT2)

!  Taylor series for computing the logarithm of a rational number
!  near 1.

      IMPLICIT NONE

!             Scratch array usage during FMLNI2:   M01 - M02

! .. Intrinsic Functions ..
      INTRINSIC int, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: int1, int2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: j, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmdivi, fmeq, fmi2m, fmmpyi
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      CALL fmi2m(int1,m02)
      CALL fmdivi(m02,int2,m02)
      CALL fmeq(m02,ma)
      ndsave = ndig
      j = 1

10    j = j + 1
      IF (int1/=1) CALL fmmpyi(m02,int1,m02)
      CALL fmdivi(m02,int2,m02)
      CALL fmdivi(m02,j,m01)
      ndig = ndsave
      CALL fmadd(ma,m01,ma)
      ndig = ndsave - int(ma(1)-m01(1))
      IF (ndig<2) ndig = 2
      IF (kflag/=1) GO TO 10

      ndig = ndsave
      ma(0) = nint(ndig*alogm2)
      ma(2) = -ma(2)
      RETURN
    END SUBROUTINE fmlni2
    SUBROUTINE fmm2dp(ma,x)

!  X = MA

!  Convert an FM number to double precision.

!  If KFLAG = -4 is returned for a value of MA that is in the range
!  of the machine's double precision number system, change the
!  definition of DPMAX in routine FMSET to reflect the current machine's
!  range.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dble
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL (KIND(0.0D0)) :: x
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: kreslt
! ..
! .. External Subroutines ..
      EXTERNAL fmargs, fmmd, fmntr, fmntrr, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      namest(ncall) = 'FMM2DP'
      kreslt = 0
      IF (abs(ma(1))>mexpab) THEN
        CALL fmargs('FMM2DP',1,ma,ma,kreslt)
      END IF
      IF (ntrace/=0) CALL fmntr(2,ma,ma,1)
      IF (kreslt/=0) THEN

!             Here no valid result can be returned.  Set X to some
!             value that the user is likely to recognize as wrong.

        x = dble(runkno)
        kflag = -4
        IF (ma(1)/=munkno) CALL fmwarn
        IF (ntrace/=0) CALL fmntrr(1,x,1)
        ncall = ncall - 1
        RETURN
      END IF

      CALL fmmd(ma,x)

      IF (ntrace/=0) CALL fmntrr(1,x,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmm2dp
    SUBROUTINE fmm2i(ma,ival)

!  IVAL = MA

!  Convert an FM number to integer.

!  KFLAG =  0 is returned if the conversion is exact.
!        = -4 is returned if MA is larger than INTMAX in magnitude.
!             IVAL = IUNKNO is returned as an indication that IVAL
!             could not be computed without integer overflow.
!        =  2 is returned if MA is smaller than INTMAX in magnitude
!             but MA is not an integer.  The next integer toward zero
!             is returned in IVAL.
!  It is sometimes convenient to call FMM2I to see if an FM number
!  can be represented as a one-word integer, by checking KFLAG upon
!  return.  To avoid an unwanted error message being printed in the
!  KFLAG=-4 case, set KWARN=0 before the call to FMM2I and reset it
!  after the call.

!  This routine performs the trace printing for the conversion.
!  FMMI is used to do the arithmetic.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. External Subroutines ..
      EXTERNAL fmmi, fmntr, fmntri
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      namest(ncall) = 'FMM2I '
      IF (ntrace/=0) CALL fmntr(2,ma,ma,1)

      CALL fmmi(ma,ival)

      IF (ntrace/=0) CALL fmntri(1,ival,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmm2i
    SUBROUTINE fmm2sp(ma,x)

!  X = MA

!  Convert an FM number to single precision.

!  MA is converted and the result is returned in X.

!  If KFLAG = -4 is returned for a value of MA that is in the range
!  of the machine's single precision number system, change the
!  definition of SPMAX in routine FMSET to reflect the current machine's
!  range.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dble, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL :: x
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: y
      INTEGER :: kreslt
! ..
! .. External Subroutines ..
      EXTERNAL fmargs, fmmd, fmntr, fmntrr, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      namest(ncall) = 'FMM2SP'
      kreslt = 0
      IF (abs(ma(1))>mexpab) THEN
        CALL fmargs('FMM2SP',1,ma,ma,kreslt)
      END IF
      IF (ntrace/=0) CALL fmntr(2,ma,ma,1)
      IF (kreslt/=0) THEN

!             Here no valid result can be returned.  Set X to some
!             value that the user is likely to recognize as wrong.

        x = runkno
        kflag = -4
        IF (ma(1)/=munkno) CALL fmwarn
        y = dble(x)
        IF (ntrace/=0) CALL fmntrr(1,y,1)
        ncall = ncall - 1
        RETURN
      END IF

      CALL fmmd(ma,y)
      x = real(y)

      IF (ntrace/=0) THEN
        y = dble(x)
        CALL fmntrr(1,y,1)
      END IF
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmm2sp
    SUBROUTINE fmmax(ma,mb,mc)

!  MC = MAX(MA,MB)

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: kwrnsv
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL fmeq, fmim, fmntr
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      kflag = 0
      ncall = ncall + 1
      namest(ncall) = 'FMMAX '
      IF (ntrace/=0) CALL fmntr(2,ma,mb,2)

      kwrnsv = kwarn
      kwarn = 0
      IF (ma(1)==munkno .OR. mb(1)==munkno) THEN
        CALL fmim(0,mc)
        mc(1) = munkno
        mc(2) = 1
        mc(0) = nint(ndig*alogm2)
        kflag = -4
      ELSE IF (fmcomp(ma,'LT',mb)) THEN
        CALL fmeq(mb,mc)
      ELSE
        CALL fmeq(ma,mc)
      END IF

      kwarn = kwrnsv
      IF (ntrace/=0) CALL fmntr(1,mc,mc,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmmax
    SUBROUTINE fmmd(ma,x)

!  X = MA

!  Internal routine for conversion to double precision.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dble, log
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL (KIND(0.0D0)) :: x
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma1, ma2
      REAL (KIND(0.0D0)) :: dlogdp, one, pmax, rzero, xbase, y, yt
      INTEGER :: j, kwrnsv, n1, ncase
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmmi, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
!             Check to see if MA is in range for single or double
!             precision.

      IF (mblogs/=mbase) CALL fmcons
      pmax = dpmax
      IF (ncall>0) THEN
        IF (namest(ncall)=='FMM2SP') pmax = dble(spmax)
      END IF
      dlogdp = log(pmax)
      ma1 = ma(1)
      ncase = 0
      IF (dble(ma(1)-1)*dlogmb>dlogdp) THEN
        kflag = -4
        x = dble(runkno)
        CALL fmwarn
        RETURN
      ELSE IF (dble(ma(1)+1)*dlogmb>dlogdp) THEN
        ma(1) = ma(1) - 2
        ncase = 1
      ELSE IF (dble(ma(1)+1)*dlogmb<-dlogdp) THEN
        kflag = -10
        x = 0.0D0
        CALL fmwarn
        RETURN
      ELSE IF (dble(ma(1)-1)*dlogmb<-dlogdp) THEN
        ma(1) = ma(1) + 2
        ncase = 2
      END IF

!             Try FMMI first so that small integers will be
!             converted exactly.

      kwrnsv = kwarn
      kwarn = 0
      CALL fmmi(ma,j)
      kwarn = kwrnsv
      IF (kflag==0) THEN
        x = j
        RETURN
      END IF
      kflag = 0

      ma2 = ma(2)
      ma(2) = abs(ma2)
      rzero = 0.0D0
      one = 1.0D0
      n1 = ndig + 1
      xbase = mbase
      x = rzero
      y = one
      DO 10 j = 2, n1
        y = y/xbase
        yt = ma(j)
        x = x + y*yt
        yt = one + y*xbase
        IF (yt<=one) GO TO 20
10    CONTINUE

20    x = x*xbase**ma(1)
      IF (ma2<0) x = -x
      ma(2) = ma2

!             Check the result if it is near overflow or underflow.

      IF (ncase==1) THEN
        IF (x<=pmax/(xbase*xbase)) THEN
          x = x*xbase*xbase
        ELSE
          kflag = -4
          x = dble(runkno)
          CALL fmwarn
        END IF
      ELSE IF (ncase==2) THEN
        IF (x>=(1.0D0/pmax)*xbase*xbase) THEN
          x = x/(xbase*xbase)
        ELSE
          kflag = -10
          x = 0.0D0
          CALL fmwarn
        END IF
      END IF
      ma(1) = ma1
      RETURN
    END SUBROUTINE fmmd
    SUBROUTINE fmmi(ma,ival)

!  IVAL = MA.  Internal FM to integer conversion routine.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, int
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: j, ka, kb, large, n1
! ..
! .. External Subroutines ..
      EXTERNAL fmwarn
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      kflag = 0
      n1 = ndig + 1
      large = int(intmax/mbase)
      ival = 0
      IF (ma(1)<=0) THEN
        IF (ma(2)/=0) kflag = 2
        RETURN
      END IF

      kb = int(ma(1)) + 1
      ival = int(abs(ma(2)))
      IF (kb>=3) THEN
        DO 10 j = 3, kb
          IF (ival>large) THEN
            kflag = -4
            IF (ma(1)/=munkno) CALL fmwarn
            ival = iunkno
            RETURN
          END IF
          IF (j<=n1) THEN
            ival = ival*int(mbase)
            IF (ival>intmax-ma(j)) THEN
              kflag = -4
              IF (ma(1)/=munkno) CALL fmwarn
              ival = iunkno
              RETURN
            ELSE
              ival = ival + int(ma(j))
            END IF
          ELSE
            ival = ival*int(mbase)
          END IF
10      CONTINUE
      END IF

      IF (ma(2)<0) ival = -ival

!             Check to see if MA is an integer.

      ka = kb + 1
      IF (ka<=n1) THEN
        DO 20 j = ka, n1
          IF (ma(j)/=0) THEN
            kflag = 2
            RETURN
          END IF
20      CONTINUE
      END IF

      RETURN
    END SUBROUTINE fmmi
    SUBROUTINE fmmin(ma,mb,mc)

!  MC = MIN(MA,MB)

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: kwrnsv
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL fmeq, fmim, fmntr
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      kflag = 0
      ncall = ncall + 1
      namest(ncall) = 'FMMIN '
      IF (ntrace/=0) CALL fmntr(2,ma,mb,2)

      kwrnsv = kwarn
      kwarn = 0
      IF (ma(1)==munkno .OR. mb(1)==munkno) THEN
        CALL fmim(0,mc)
        mc(1) = munkno
        mc(2) = 1
        mc(0) = nint(ndig*alogm2)
        kflag = -4
      ELSE IF (fmcomp(ma,'GT',mb)) THEN
        CALL fmeq(mb,mc)
      ELSE
        CALL fmeq(ma,mc)
      END IF

      kwarn = kwrnsv
      IF (ntrace/=0) CALL fmntr(1,mc,mc,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmmin
    SUBROUTINE fmmod(ma,mb,mc)

!  MC = MA(MOD MB).

      IMPLICIT NONE

!             Scratch array usage during FMMOD:   M01 - M03

! .. Intrinsic Functions ..
      INTRINSIC abs, int, log, max, min, mod, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: macca, maccb, macmax, mvb, mvc, mvy, mvz, mxsave
      INTEGER :: j, k, kasave, kb, ke, kn, kovun, kreslt, kwrnsv, ndsave, &
        ntrsav
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmdiv, fmentr, fmeq2, fmexit, fmi2m, fmint, &
        fmm2i, fmmpy, fmntr, fmrslt, fmsub, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (abs(ma(1))>mexpab .OR. abs(mb(1))>mexpab) THEN
        CALL fmentr('FMMOD ',ma,mb,2,mc,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        namest(ncall) = 'FMMOD '
        IF (ntrace/=0) CALL fmntr(2,ma,mb,2)
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
        IF (mb(1)==mexpov .OR. mb(1)==mexpun) kovun = 1
        ndsave = ndig
        IF (ncall==1) THEN
          k = max(ngrd52-1,2)
          ndig = max(ndig+k,2)
          IF (ndig>ndg2mx) THEN
            kflag = -9
            CALL fmwarn
            ndig = ndsave
            kreslt = 12
            CALL fmrslt(ma,mb,mc,kreslt)
            IF (ntrace/=0) CALL fmntr(1,mc,mc,1)
            ncall = ncall - 1
            RETURN
          END IF
        END IF
        kasave = kaccsw
        mxsave = mxexp
        mxexp = mxexp2
      END IF
      kwrnsv = kwarn
      kwarn = 0
      macca = ma(0)
      maccb = mb(0)

      IF (mb(1)>ma(1) .AND. mb(2)/=0) THEN
        CALL fmeq2(ma,m01,ndsave,ndig,0)
        m01(0) = nint(ndig*alogm2)
      ELSE

!             Special cases when MB is a small integer.

        CALL fmeq2(ma,m02,ndsave,ndig,0)
        m02(0) = nint(ndig*alogm2)
        CALL fmeq2(mb,m03,ndsave,ndig,0)
        m03(0) = nint(ndig*alogm2)

        CALL fmm2i(m03,kb)
        IF (kflag==0 .AND. kb<mxbase) THEN
          IF (kb==1 .OR. kb==-1) THEN
            IF (m02(1)>=ndig) THEN
              CALL fmi2m(0,m01)
              GO TO 70
            ELSE
              CALL fmint(m02,m03)
              CALL fmsub(m02,m03,m01)
              GO TO 70
            END IF
          ELSE IF (m02(1)==mexpov .OR. kb==0) THEN
            kflag = -4
            kwarn = kwrnsv
            kaccsw = kasave
            mxexp = mxsave
            CALL fmwarn
            mc(1) = munkno
            mc(2) = 1
            mc(0) = nint(ndig*alogm2)
            DO 10 j = 2, ndsave
              mc(j+1) = 0
10          CONTINUE
            ndig = ndsave
            IF (ntrace/=0) CALL fmntr(1,mc,mc,1)
            ncall = ncall - 1
            RETURN
          ELSE IF (m02(1)>ndig .AND. mod(int(mbase),kb)==0) THEN
            CALL fmi2m(0,m01)
            GO TO 70
          END IF
          IF (m02(1)<ndig) THEN
            DO 20 j = int(m02(1)) + 1, ndig + 1
              IF (m02(j)/=0) GO TO 50
20          CONTINUE
          END IF
          ke = min(int(m02(1)),ndig)
          mvb = kb
          mvc = mod(m02(2),mvb)
          DO 30 j = 3, ke + 1
            mvc = mod(mvc*mbase+m02(j),mvb)
30        CONTINUE
          IF (mvc==0) THEN
            CALL fmi2m(0,m01)
            GO TO 70
          END IF
          kn = int(m02(1)) - ke
          mvy = mod(mbase,mvb)
          mvz = 1
          IF (mod(kn,2)==1) mvz = mvy

          IF (mvy/=1) THEN
40          kn = kn/2
            mvy = mod(mvy*mvy,mvb)
            IF (mod(kn,2)==1) mvz = mod(mvz*mvy,mvb)
            IF (kn>1) GO TO 40
          END IF
          mvz = mod(mvz*mvc,mvb)
          ke = int(mvz)
          CALL fmi2m(ke,m01)
          GO TO 70
        END IF

!             General case.

50      IF (ma(2)/=0) THEN
          ndig = ndig + int(ma(1)-mb(1))
        END IF
        IF (ndig>ndg2mx .OR. mb(2)==0) THEN
          kflag = -9
          IF (ma(1)==mexpov .OR. mb(1)==mexpun .OR. mb(2)==0) kflag = -4
          kwarn = kwrnsv
          kaccsw = kasave
          mxexp = mxsave
          CALL fmwarn
          mc(1) = munkno
          mc(2) = 1
          mc(0) = nint(ndig*alogm2)
          DO 60 j = 2, ndsave
            mc(j+1) = 0
60        CONTINUE
          ndig = ndsave
          IF (ntrace/=0) CALL fmntr(1,mc,mc,1)
          ncall = ncall - 1
          RETURN
        END IF

        CALL fmeq2(ma,m02,ndsave,ndig,0)
        m02(0) = nint(ndig*alogm2)
        CALL fmeq2(mb,m03,ndsave,ndig,0)
        m03(0) = nint(ndig*alogm2)

        m02(2) = abs(m02(2))
        m03(2) = abs(m03(2))
        CALL fmdiv(m02,m03,m01)
        CALL fmint(m01,m01)
        CALL fmmpy(m01,m03,m01)
        CALL fmsub(m02,m01,m01)

!             Due to rounding, M01 may not be between 0 and MB here.

        ntrsav = ntrace
        ntrace = 0
        IF (fmcomp(m01,'GE',m03)) THEN
          ntrace = ntrsav
          CALL fmsub(m01,m03,m01)
        END IF
        ntrace = ntrsav
        IF (m01(2)<0) CALL fmadd(m01,m03,m01)
        IF (ma(2)<0 .AND. m01(1)/=munkno) m01(2) = -m01(2)
      END IF

70    IF (kflag==1) kflag = 0
      kwarn = kwrnsv
      macmax = nint((ndsave-1)*alogm2+log(real(abs(m01(2))+1))/0.69315)
      m01(0) = min(m01(0),macca,maccb,macmax)
      CALL fmexit(m01,mc,ndsave,mxsave,kasave,kovun)
      RETURN
    END SUBROUTINE fmmod
    SUBROUTINE fmmove(mw,ma)

!  Move a result from a work area (MW) to MA.

!  If the result has MW(2)=0, then it is shifted and the exponent
!  adjusted when it is moved to MA.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mw(lmwa)
! ..
! .. Local Scalars ..
      INTEGER :: j, n1, n2
! ..
! .. External Subroutines ..
      EXTERNAL fmtrap
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mw(2)/=0) THEN
        n1 = ndig + 1

!             Major (Inner Loop)

        DO 10 j = 1, n1
          ma(j) = mw(j)
10      CONTINUE
      ELSE
        n2 = ndig + 2
        DO 20 j = 3, n2
          ma(j-1) = mw(j)
20      CONTINUE
        IF (ma(2)/=0) THEN
          ma(1) = mw(1) - 1
        ELSE
          ma(1) = 0
        END IF
      END IF

      IF (abs(ma(1))>mxexp) CALL fmtrap(ma)

      RETURN
    END SUBROUTINE fmmove
    SUBROUTINE fmmpy(ma,mb,mc)

!  MC = MA * MB

!  When one of the numbers MA, MB is known to have more zero digits
!  (base MBASE) than the other, it is faster if MB is the one with
!  more zero digits.

!  This routine performs the trace printing for multiplication.
!  FMMPY2 is used to do the arithmetic.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. External Subroutines ..
      EXTERNAL fmmpy2, fmntr
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (ntrace/=0) THEN
        namest(ncall) = 'FMMPY '
        CALL fmntr(2,ma,mb,2)

        CALL fmmpy2(ma,mb,mc)

        CALL fmntr(1,mc,mc,1)
      ELSE
        CALL fmmpy2(ma,mb,mc)
      END IF
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmmpy
    SUBROUTINE fmmpy2(ma,mb,mc)

!  Internal multiplication routine.  MC = MA * MB

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, log, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, macca, maccb, mb2, md2b, mr
      INTEGER :: j, kreslt, kshift, n1, nguard, nzma, nzmb
! ..
! .. External Subroutines ..
      EXTERNAL fmargs, fmcons, fmim, fmmove, fmmpy3, fmrnd, fmrslt, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      macca = ma(0)
      maccb = mb(0)
      IF (abs(ma(1))>mexpab .OR. abs(mb(1))>mexpab .OR. kdebug==1) THEN
        CALL fmargs('FMMPY ',2,ma,mb,kreslt)
        IF (kreslt/=0) THEN
          ncall = ncall + 1
          namest(ncall) = 'FMMPY '
          CALL fmrslt(ma,mb,mc,kreslt)
          ncall = ncall - 1
          RETURN
        END IF
      ELSE IF (ma(2)==0 .OR. mb(2)==0) THEN
        CALL fmim(0,mc)
        mc(0) = min(macca,maccb)
        RETURN
      END IF
      kflag = 0

!             Save the sign of MA and MB and then work only with
!             positive numbers.

      ma2 = ma(2)
      mb2 = mb(2)
      ma(2) = abs(ma(2))
      mb(2) = abs(mb(2))

!             NGUARD is the number of guard digits used.

      IF (ncall>1) THEN
        nguard = ngrd22
        IF (nguard>ndig) nguard = ndig
      ELSE
        nguard = ngrd52
        IF (nguard>ndig) nguard = ndig
      END IF
      IF (ma(2)*mb(2)<mbase .AND. nguard<3) nguard = 3

      n1 = ndig + 1

!             If there is a good chance of finding several zero digits,
!             see which number has more zero digits.

      IF (ndig>=6*mbase) THEN
        nzma = 0
        nzmb = 0
        DO 10 j = 2, n1
          IF (ma(j)==0) nzma = nzma + 1
          IF (mb(j)==0) nzmb = nzmb + 1
10      CONTINUE

!             It is faster if the second argument is the one with
!             more zero digits.

        IF (nzma>nzmb) THEN
          CALL fmmpy3(mb,ma,nguard,kshift)
        ELSE
          CALL fmmpy3(ma,mb,nguard,kshift)
        END IF
      ELSE
        CALL fmmpy3(ma,mb,nguard,kshift)
      END IF

!             The multiplication is complete.  Round the result,
!             move it to MC, and append the correct sign.

      ma(2) = ma2
      mb(2) = mb2
      mr = 2*mwa(ndig+2+kshift) + 1
      IF (mr>=mbase) THEN
        IF (mr-1>mbase .AND. mwa(n1+kshift)<mbase-1) THEN
          IF (kround/=0 .OR. ncall>1) THEN
            mwa(n1+kshift) = mwa(n1+kshift) + 1
            mwa(n1+1+kshift) = 0
          END IF
        ELSE
          CALL fmrnd(mwa,ndig,nguard,kshift)
        END IF
      END IF
      CALL fmmove(mwa,mc)

      IF (kflag<0) THEN
        namest(ncall) = 'FMMPY '
        CALL fmwarn
      END IF

      IF (ma2*mb2<0) mc(2) = -mc(2)

      IF (kaccsw==1) THEN
        md2b = nint((ndig-1)*alogm2+log(real(abs(mc(2))+1))/0.69315)
        mc(0) = min(macca,maccb,md2b)
      ELSE
        mc(0) = min(macca,maccb)
      END IF
      RETURN
    END SUBROUTINE fmmpy2
    SUBROUTINE fmmpy3(ma,mb,nguard,kshift)

!  Internal multiplication of MA*MB.  The result is returned in MWA.
!  Both MA and MB are positive.

!  NGUARD is the number of guard digits that will be used.
!  KSHIFT = 1 is returned if a left shift is pending (i.e., MWA(2)=0).
!             The shift will be done in FMMOVE.  KSHIFT = 0 is returned
!             if no shift is pending.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC dint, int, min
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kshift, nguard
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maxmwa, mbj, mbkj, mbm1, mbnorm, mk, mkt, mmax, mt
      INTEGER :: j, jm1, k, kb, ki, kj, kl, knz, kwa, l, n1
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      n1 = ndig + 1
      mwa(1) = ma(1) + mb(1)
      l = n1 + nguard
      mwa(l+1) = 0

!             The multiplication loop begins here.

!             MBNORM is the minimum number of digits that can be
!                    multiplied before normalization is required.
!             MAXMWA is an upper bound on the size of values in MWA
!                    divided by (MBASE-1).  It is used to determine
!                    whether to normalize before the next digit is
!                    multiplied.

      mbm1 = mbase - 1
      mbnorm = dint(maxint/(mbm1*mbm1))
      mmax = intmax - mbase
      mmax = min(dint(maxint/mbm1-mbm1),mmax)
      IF (mbnorm>1) THEN
        mbj = mb(2)

!             Count the trailing zeros in MA.

        IF (ma(n1)/=0) THEN
          knz = n1
        ELSE
          DO 10 j = ndig, 2, -1
            IF (ma(j)/=0) THEN
              knz = j
              GO TO 20
            END IF
10        CONTINUE
        END IF

20      mwa(2) = 0
        DO 30 k = ndig + 2, l
          mwa(k) = 0
30      CONTINUE

!             (Inner Loop)

        DO 40 k = 2, n1
          mwa(k+1) = ma(k)*mbj
40      CONTINUE
        maxmwa = mbj
        DO 70 j = 3, n1
          mbj = mb(j)
          IF (mbj/=0) THEN
            maxmwa = maxmwa + mbj
            jm1 = j - 1
            kl = min(knz,l-jm1)

!                       Major (Inner Loop)

            DO 50 k = j + 1, j + kl - 1
              mwa(k) = mwa(k) + ma(k-jm1)*mbj
50          CONTINUE
          END IF

          IF (maxmwa>mmax) THEN
            maxmwa = 0

!                       Here normalization is only required for the
!                       range of digits currently changing in MWA.

            DO 60 kb = jm1 + kl, jm1 + 2, -1
              mkt = int(mwa(kb)/mbase)
              mwa(kb-1) = mwa(kb-1) + mkt
              mwa(kb) = mwa(kb) - mkt*mbase
60          CONTINUE
          END IF
70      CONTINUE

!             Perform the final normalization.  (Inner Loop)

        DO 80 kb = l, 3, -1
          mkt = int(mwa(kb)/mbase)
          mwa(kb-1) = mwa(kb-1) + mkt
          mwa(kb) = mwa(kb) - mkt*mbase
80      CONTINUE

      ELSE

!             If normalization must be done for each digit, combine
!             the two loops and normalize as the digits are multiplied.

        DO 90 j = 2, l
          mwa(j) = 0
90      CONTINUE
        kj = ndig + 2
        DO 110 j = 2, n1
          kj = kj - 1
          mbkj = mb(kj)
          IF (mbkj==0) GO TO 110
          kl = l - kj + 1
          IF (kl>n1) kl = n1
          ki = kl + 2
          kwa = kl + kj + 1
          mk = 0
          DO 100 k = 2, kl
            mt = ma(ki-k)*mbkj + mwa(kwa-k) + mk
            mk = int(mt/mbase)
            mwa(kwa-k) = mt - mbase*mk
100       CONTINUE
          mwa(kwa-kl-1) = mk
110     CONTINUE

      END IF

!             Set KSHIFT = 1 if a shift left is necessary.

      IF (mwa(2)==0) THEN
        kshift = 1
        RETURN
      ELSE
        kshift = 0
        RETURN
      END IF

    END SUBROUTINE fmmpy3
    SUBROUTINE fmmpyd(ma,mb,mc,md,me)

!  Double multiplication routine.  MD = MA * MB,   ME = MA * MC

!  It is usually slightly faster to do two multiplications that
!  have a common factor with one call.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dint, int, log, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck), &
        me(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, macca, maccb, maccc, maxmwa, mb2, mbj, mbkj, mbm1, &
        mbnorm, mc2, mcj, mckj, md2b, mkb, mkc, mkt, mmax, mr, mt, mtemp
      INTEGER :: j, jm1, k, kb, ki, kj, kl, knz, kovun, kshift, kwa, l, n1, &
        nguard
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmeq, fmim, fmmove, fmmpy2, fmntr, fmntrj, fmprnt, &
        fmrnd, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa), mwd(lmwa), mwe(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
      COMMON /fmwa/mwd, mwe
! ..
      ncall = ncall + 1
      IF (ntrace/=0) THEN
        namest(ncall) = 'FMMPYD'
        CALL fmntr(2,ma,mb,2)
        IF (abs(ntrace)>=2 .AND. ncall<=lvltrc) THEN
          IF (ntrace<0) THEN
            CALL fmntrj(mc,ndig)
          ELSE
            CALL fmprnt(mc)
          END IF
        END IF
      END IF

      IF (mblogs/=mbase) CALL fmcons
      macca = ma(0)
      maccb = mb(0)
      maccc = mc(0)
      IF (abs(ma(1))>mexpab .OR. abs(mb(1))>mexpab .OR. abs(mc(1))>mexpab) &
          THEN
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun .OR. mb(1)==mexpov .OR. &
          mb(1)==mexpun .OR. mc(1)==mexpov .OR. mc(1)==mexpun) kovun = 1
        IF (ma(1)==munkno .OR. mb(1)==munkno .OR. mc(1)==munkno) kovun = 2
        ncall = ncall + 1
        CALL fmmpy2(ma,mb,mwd)
        kb = kflag
        CALL fmmpy2(ma,mc,me)
        ncall = ncall - 1
        IF (((kflag<0 .OR. kb<0) .AND. kovun==0) .OR. ((kflag==-4 .OR. kb== &
            -4) .AND. kovun==1)) THEN
          IF (kflag==-4 .OR. kb==-4) THEN
            kflag = -4
          ELSE IF (kflag==-5 .OR. kb==-5) THEN
            kflag = -5
          ELSE
            kflag = min(kflag,kb)
          END IF
          namest(ncall) = 'FMMPYD'
          CALL fmwarn
        END IF
        CALL fmeq(mwd,md)
        GO TO 120
      END IF
      IF (ma(2)==0) THEN
        CALL fmim(0,md)
        md(0) = min(macca,maccb)
        CALL fmim(0,me)
        me(0) = min(macca,maccc)
        GO TO 120
      END IF
      IF (mb(2)==0) THEN
        CALL fmmpy2(ma,mc,me)
        CALL fmim(0,md)
        md(0) = min(macca,maccb)
        GO TO 120
      END IF
      IF (mc(2)==0) THEN
        CALL fmmpy2(ma,mb,md)
        CALL fmim(0,me)
        me(0) = min(macca,maccc)
        GO TO 120
      END IF
      kflag = 0

!             NGUARD is the number of guard digits used.

      IF (ncall>1) THEN
        nguard = ngrd22
        IF (nguard>ndig) nguard = ndig
      ELSE
        nguard = ngrd52
        IF (nguard>ndig) nguard = ndig
      END IF
      IF ((ma(2)*mb(2)<mbase .OR. ma(2)*mc(2)<mbase) .AND. nguard<3) &
        nguard = 3

!             Save the sign of MA, MB, and MC and then
!             work only with positive numbers.

      ma2 = ma(2)
      mb2 = mb(2)
      mc2 = mc(2)
      ma(2) = abs(ma(2))
      mb(2) = abs(mb(2))
      mc(2) = abs(mc(2))

      n1 = ndig + 1
      mwa(1) = ma(1) + mb(1)
      mwd(1) = ma(1) + mc(1)
      l = ndig + 1 + nguard
      mwa(l+1) = 0
      mwd(l+1) = 0

!             The multiplication loop begins here.

!             MBNORM is the minimum number of digits that can be
!                    multiplied before normalization is required.
!             MAXMWA is an upper bound on the size of values in MWA
!                    divided by (MBASE-1).  It is used to determine
!                    whether to normalize before the next digit is
!                    multiplied.

      mbm1 = mbase - 1
      mbnorm = dint(maxint/(mbm1*mbm1))
      mmax = intmax - mbase
      mmax = min(dint(maxint/mbm1-mbm1),mmax)
      IF (mbnorm>1) THEN
        mbj = mb(2)
        mcj = mc(2)

!             Count the trailing zeros in MA.

        IF (ma(n1)/=0) THEN
          knz = n1
        ELSE
          DO 10 j = ndig, 2, -1
            IF (ma(j)/=0) THEN
              knz = j
              GO TO 20
            END IF
10        CONTINUE
        END IF

20      mwa(2) = 0
        mwd(2) = 0
        DO 30 k = ndig + 2, l
          mwa(k) = 0
          mwd(k) = 0
30      CONTINUE

!             (Inner Loop)

        DO 40 k = 2, n1
          mtemp = ma(k)
          mwa(k+1) = mtemp*mbj
          mwd(k+1) = mtemp*mcj
40      CONTINUE
        IF (mbj>mcj) THEN
          maxmwa = mbj
        ELSE
          maxmwa = mcj
        END IF
        DO 70 j = 3, n1
          mbj = mb(j)
          mcj = mc(j)
          IF (mbj>mcj) THEN
            maxmwa = maxmwa + mbj
          ELSE
            maxmwa = maxmwa + mcj
          END IF
          jm1 = j - 1
          kl = min(knz,l-jm1)

!                       Major (Inner Loop)

          DO 50 k = j + 1, j + kl - 1
            mtemp = ma(k-jm1)
            mwa(k) = mwa(k) + mtemp*mbj
            mwd(k) = mwd(k) + mtemp*mcj
50        CONTINUE

          IF (maxmwa>mmax) THEN
            maxmwa = 0

!                       Here normalization is only required for the
!                       range of digits currently changing in MWA.

            DO 60 kb = jm1 + kl, jm1 + 2, -1
              mkt = int(mwa(kb)/mbase)
              mwa(kb-1) = mwa(kb-1) + mkt
              mwa(kb) = mwa(kb) - mkt*mbase
              mkt = int(mwd(kb)/mbase)
              mwd(kb-1) = mwd(kb-1) + mkt
              mwd(kb) = mwd(kb) - mkt*mbase
60          CONTINUE
          END IF
70      CONTINUE

!             Perform the final normalization.  (Inner Loop)

        DO 80 kb = l, 3, -1
          mkt = int(mwa(kb)/mbase)
          mwa(kb-1) = mwa(kb-1) + mkt
          mwa(kb) = mwa(kb) - mkt*mbase
          mkt = int(mwd(kb)/mbase)
          mwd(kb-1) = mwd(kb-1) + mkt
          mwd(kb) = mwd(kb) - mkt*mbase
80      CONTINUE

      ELSE

!             If normalization must be done for each digit, combine
!             the two loops and normalize as the digits are multiplied.

        DO 90 j = 2, l
          mwa(j) = 0
          mwd(j) = 0
90      CONTINUE
        kj = ndig + 2
        DO 110 j = 2, n1
          kj = kj - 1
          mbkj = mb(kj)
          mckj = mc(kj)
          kl = l - kj + 1
          IF (kl>n1) kl = n1
          ki = kl + 2
          kwa = kl + kj + 1
          mkb = 0
          mkc = 0
          DO 100 k = 2, kl
            mt = ma(ki-k)*mbkj + mwa(kwa-k) + mkb
            mkb = int(mt/mbase)
            mwa(kwa-k) = mt - mbase*mkb
            mt = ma(ki-k)*mckj + mwd(kwa-k) + mkc
            mkc = int(mt/mbase)
            mwd(kwa-k) = mt - mbase*mkc
100       CONTINUE
          mwa(kwa-kl-1) = mkb
          mwd(kwa-kl-1) = mkc
110     CONTINUE

      END IF

!             Set KSHIFT = 1 if a shift left is necessary.

      IF (mwa(2)==0) THEN
        kshift = 1
      ELSE
        kshift = 0
      END IF

!             The multiplications are complete.

      ma(2) = ma2
      mb(2) = mb2
      mc(2) = mc2
      mr = 2*mwa(ndig+2+kshift) + 1
      IF (mr>=mbase) THEN
        IF (mr-1>mbase .AND. mwa(n1+kshift)<mbase-1) THEN
          IF (kround/=0 .OR. ncall>1) THEN
            mwa(n1+kshift) = mwa(n1+kshift) + 1
            mwa(n1+1+kshift) = 0
          END IF
        ELSE
          CALL fmrnd(mwa,ndig,nguard,kshift)
        END IF
      END IF
      CALL fmmove(mwa,md)

      IF (mwd(2)==0) THEN
        kshift = 1
      ELSE
        kshift = 0
      END IF
      mr = 2*mwd(ndig+2+kshift) + 1
      IF (mr>=mbase) THEN
        IF (mr-1>mbase .AND. mwd(n1+kshift)<mbase-1) THEN
          IF (kround/=0 .OR. ncall>1) THEN
            mwd(n1+kshift) = mwd(n1+kshift) + 1
            mwd(n1+1+kshift) = 0
          END IF
        ELSE
          CALL fmrnd(mwd,ndig,nguard,kshift)
        END IF
      END IF
      CALL fmmove(mwd,me)

      IF (kflag<0) THEN
        namest(ncall) = 'FMMPYD'
        CALL fmwarn
      END IF

      IF (ma2*mb2<0) md(2) = -md(2)
      IF (ma2*mc2<0) me(2) = -me(2)

      IF (kaccsw==1) THEN
        md2b = nint((ndig-1)*alogm2+log(real(abs(md(2))+1))/0.69315)
        md(0) = min(macca,maccb,md2b)
        md2b = nint((ndig-1)*alogm2+log(real(abs(me(2))+1))/0.69315)
        me(0) = min(macca,maccc,md2b)
      ELSE
        md(0) = min(macca,maccb)
        me(0) = min(macca,maccc)
      END IF

120   IF (ntrace/=0) THEN
        CALL fmntr(1,md,md,1)
        IF (abs(ntrace)>=1 .AND. ncall<=lvltrc) THEN
          IF (ntrace<0) THEN
            CALL fmntrj(me,ndig)
          ELSE
            CALL fmprnt(me)
          END IF
        END IF
      END IF
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmmpyd
    SUBROUTINE fmmpye(ma,mb,mc,md,me,mf,mg)

!  Triple multiplication routine.

!      ME = MA * MB,   MF = MA * MC,   MG = MA * MD

!  It is usually slightly faster to do three multiplications that
!  have a common factor with one call.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dint, int, log, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck), &
        me(0:lunpck), mf(0:lunpck), mg(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, macca, maccb, maccc, maccd, maxj, maxmwa, mb2, mbj, &
        mbkj, mbm1, mbnorm, mc2, mcj, mckj, md2, md2b, mdj, mdkj, mkb, mkc, &
        mkd, mkt, mmax, mr, mt, mtemp
      INTEGER :: j, jm1, k, kb, ki, kj, kl, knz, kovun, kshift, kwa, l, n1, &
        nguard
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmeq, fmim, fmmove, fmmpy2, fmntr, fmntrj, fmprnt, &
        fmrnd, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa), mwd(lmwa), mwe(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
      COMMON /fmwa/mwd, mwe
! ..
      ncall = ncall + 1
      IF (ntrace/=0) THEN
        namest(ncall) = 'FMMPYE'
        CALL fmntr(2,ma,mb,2)
        IF (abs(ntrace)>=2 .AND. ncall<=lvltrc) THEN
          IF (ntrace<0) THEN
            CALL fmntrj(mc,ndig)
            CALL fmntrj(md,ndig)
          ELSE
            CALL fmprnt(mc)
            CALL fmprnt(md)
          END IF
        END IF
      END IF

      IF (mblogs/=mbase) CALL fmcons
      macca = ma(0)
      maccb = mb(0)
      maccc = mc(0)
      maccd = md(0)
      IF (abs(ma(1))>mexpab .OR. abs(mb(1))>mexpab .OR. abs(mc( &
          1))>mexpab .OR. abs(md(1))>mexpab) THEN
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun .OR. mb(1)==mexpov .OR. &
          mb(1)==mexpun .OR. mc(1)==mexpov .OR. mc(1)==mexpun .OR. &
          md(1)==mexpov .OR. md(1)==mexpun) kovun = 1
        IF (ma(1)==munkno .OR. mb(1)==munkno .OR. mc(1)==munkno .OR. &
          md(1)==munkno) kovun = 2
        ncall = ncall + 1
        CALL fmmpy2(ma,mb,mwd)
        kb = kflag
        CALL fmmpy2(ma,mc,mwe)
        kj = kflag
        CALL fmmpy2(ma,md,mg)
        ncall = ncall - 1
        IF (((kflag<0 .OR. kb<0 .OR. kj<0) .AND. kovun==0) .OR. ((kflag== &
            -4 .OR. kb==-4 .OR. kj==-4) .AND. kovun==1)) THEN
          IF (kflag==-4 .OR. kb==-4 .OR. kj==-4) THEN
            kflag = -4
          ELSE IF (kflag==-5 .OR. kb==-5 .OR. kj==-5) THEN
            kflag = -5
          ELSE
            kflag = min(kflag,kb,kj)
          END IF
          namest(ncall) = 'FMMPYE'
          CALL fmwarn
        END IF
        CALL fmeq(mwd,me)
        CALL fmeq(mwe,mf)
        GO TO 120
      END IF
      IF (ma(2)==0) THEN
        CALL fmim(0,me)
        me(0) = min(macca,maccb)
        CALL fmim(0,mf)
        mf(0) = min(macca,maccc)
        CALL fmim(0,mg)
        mg(0) = min(macca,maccd)
        GO TO 120
      END IF
      IF (mb(2)==0 .OR. mc(2)==0 .OR. md(2)==0) THEN
        CALL fmmpy2(ma,mb,mwd)
        CALL fmmpy2(ma,mc,mwe)
        CALL fmmpy2(ma,md,mg)
        CALL fmeq(mwd,me)
        CALL fmeq(mwe,mf)
        GO TO 120
      END IF
      kflag = 0

!             NGUARD is the number of guard digits used.

      IF (ncall>1) THEN
        nguard = ngrd22
        IF (nguard>ndig) nguard = ndig
      ELSE
        nguard = ngrd52
        IF (nguard>ndig) nguard = ndig
      END IF
      IF ((ma(2)*mb(2)<mbase .OR. ma(2)*mc(2)<mbase .OR. ma(2)*md( &
        2)<mbase) .AND. nguard<3) nguard = 3

!             Save the signs and then work only with positive numbers.

      ma2 = ma(2)
      mb2 = mb(2)
      mc2 = mc(2)
      md2 = md(2)
      ma(2) = abs(ma(2))
      mb(2) = abs(mb(2))
      mc(2) = abs(mc(2))
      md(2) = abs(md(2))

      n1 = ndig + 1
      mwa(1) = ma(1) + mb(1)
      mwd(1) = ma(1) + mc(1)
      mwe(1) = ma(1) + md(1)
      l = ndig + 1 + nguard
      mwa(l+1) = 0
      mwd(l+1) = 0
      mwe(l+1) = 0

!             The multiplication loop begins here.

!             MBNORM is the minimum number of digits that can be
!                    multiplied before normalization is required.
!             MAXMWA is an upper bound on the size of values in MWA
!                    divided by (MBASE-1).  It is used to determine
!                    whether to normalize before the next digit is
!                    multiplied.

      mbm1 = mbase - 1
      mbnorm = dint(maxint/(mbm1*mbm1))
      mmax = intmax - mbase
      mmax = min(dint(maxint/mbm1-mbm1),mmax)
      IF (mbnorm>1) THEN
        mbj = mb(2)
        mcj = mc(2)
        mdj = md(2)

!             Count the trailing zeros in MA.

        IF (ma(n1)/=0) THEN
          knz = n1
        ELSE
          DO 10 j = ndig, 2, -1
            IF (ma(j)/=0) THEN
              knz = j
              GO TO 20
            END IF
10        CONTINUE
        END IF

20      mwa(2) = 0
        mwd(2) = 0
        mwe(2) = 0
        DO 30 k = ndig + 2, l
          mwa(k) = 0
          mwd(k) = 0
          mwe(k) = 0
30      CONTINUE

!             (Inner Loop)

        DO 40 k = 2, n1
          mtemp = ma(k)
          mwa(k+1) = mtemp*mbj
          mwd(k+1) = mtemp*mcj
          mwe(k+1) = mtemp*mdj
40      CONTINUE
        maxmwa = mbj
        IF (mcj>maxmwa) maxmwa = mcj
        IF (mdj>maxmwa) maxmwa = mdj
        DO 70 j = 3, n1
          mbj = mb(j)
          mcj = mc(j)
          mdj = md(j)
          maxj = mbj
          IF (mcj>maxj) maxj = mcj
          IF (mdj>maxj) maxj = mdj
          maxmwa = maxmwa + maxj
          jm1 = j - 1
          kl = min(knz,l-jm1)

!                       Major (Inner Loop)

          DO 50 k = j + 1, j + kl - 1
            mtemp = ma(k-jm1)
            mwa(k) = mwa(k) + mtemp*mbj
            mwd(k) = mwd(k) + mtemp*mcj
            mwe(k) = mwe(k) + mtemp*mdj
50        CONTINUE

          IF (maxmwa>mmax) THEN
            maxmwa = 0

!                       Here normalization is only required for the
!                       range of digits currently changing in MWA.

            DO 60 kb = jm1 + kl, jm1 + 2, -1
              mkt = int(mwa(kb)/mbase)
              mwa(kb-1) = mwa(kb-1) + mkt
              mwa(kb) = mwa(kb) - mkt*mbase
              mkt = int(mwd(kb)/mbase)
              mwd(kb-1) = mwd(kb-1) + mkt
              mwd(kb) = mwd(kb) - mkt*mbase
              mkt = int(mwe(kb)/mbase)
              mwe(kb-1) = mwe(kb-1) + mkt
              mwe(kb) = mwe(kb) - mkt*mbase
60          CONTINUE
          END IF
70      CONTINUE

!             Perform the final normalization.  (Inner Loop)

        DO 80 kb = l, 3, -1
          mkt = int(mwa(kb)/mbase)
          mwa(kb-1) = mwa(kb-1) + mkt
          mwa(kb) = mwa(kb) - mkt*mbase
          mkt = int(mwd(kb)/mbase)
          mwd(kb-1) = mwd(kb-1) + mkt
          mwd(kb) = mwd(kb) - mkt*mbase
          mkt = int(mwe(kb)/mbase)
          mwe(kb-1) = mwe(kb-1) + mkt
          mwe(kb) = mwe(kb) - mkt*mbase
80      CONTINUE

      ELSE

!             If normalization must be done for each digit, combine
!             the two loops and normalize as the digits are multiplied.

        DO 90 j = 2, l
          mwa(j) = 0
          mwd(j) = 0
          mwe(j) = 0
90      CONTINUE
        kj = ndig + 2
        DO 110 j = 2, n1
          kj = kj - 1
          mbkj = mb(kj)
          mckj = mc(kj)
          mdkj = md(kj)
          kl = l - kj + 1
          IF (kl>n1) kl = n1
          ki = kl + 2
          kwa = kl + kj + 1
          mkb = 0
          mkc = 0
          mkd = 0
          DO 100 k = 2, kl
            mt = ma(ki-k)*mbkj + mwa(kwa-k) + mkb
            mkb = int(mt/mbase)
            mwa(kwa-k) = mt - mbase*mkb
            mt = ma(ki-k)*mckj + mwd(kwa-k) + mkc
            mkc = int(mt/mbase)
            mwd(kwa-k) = mt - mbase*mkc
            mt = ma(ki-k)*mdkj + mwe(kwa-k) + mkd
            mkd = int(mt/mbase)
            mwe(kwa-k) = mt - mbase*mkd
100       CONTINUE
          mwa(kwa-kl-1) = mkb
          mwd(kwa-kl-1) = mkc
          mwe(kwa-kl-1) = mkd
110     CONTINUE

      END IF

!             Set KSHIFT = 1 if a shift left is necessary.

      IF (mwa(2)==0) THEN
        kshift = 1
      ELSE
        kshift = 0
      END IF

!             The multiplications are complete.

      ma(2) = ma2
      mb(2) = mb2
      mc(2) = mc2
      md(2) = md2
      mr = 2*mwa(ndig+2+kshift) + 1
      IF (mr>=mbase) THEN
        IF (mr-1>mbase .AND. mwa(n1+kshift)<mbase-1) THEN
          IF (kround/=0 .OR. ncall>1) THEN
            mwa(n1+kshift) = mwa(n1+kshift) + 1
            mwa(n1+1+kshift) = 0
          END IF
        ELSE
          CALL fmrnd(mwa,ndig,nguard,kshift)
        END IF
      END IF
      CALL fmmove(mwa,me)

      IF (mwd(2)==0) THEN
        kshift = 1
      ELSE
        kshift = 0
      END IF
      mr = 2*mwd(ndig+2+kshift) + 1
      IF (mr>=mbase) THEN
        IF (mr-1>mbase .AND. mwd(n1+kshift)<mbase-1) THEN
          IF (kround/=0 .OR. ncall>1) THEN
            mwd(n1+kshift) = mwd(n1+kshift) + 1
            mwd(n1+1+kshift) = 0
          END IF
        ELSE
          CALL fmrnd(mwd,ndig,nguard,kshift)
        END IF
      END IF
      CALL fmmove(mwd,mf)

      IF (mwe(2)==0) THEN
        kshift = 1
      ELSE
        kshift = 0
      END IF
      mr = 2*mwe(ndig+2+kshift) + 1
      IF (mr>=mbase) THEN
        IF (mr-1>mbase .AND. mwe(n1+kshift)<mbase-1) THEN
          IF (kround/=0 .OR. ncall>1) THEN
            mwe(n1+kshift) = mwe(n1+kshift) + 1
            mwe(n1+1+kshift) = 0
          END IF
        ELSE
          CALL fmrnd(mwe,ndig,nguard,kshift)
        END IF
      END IF
      CALL fmmove(mwe,mg)

      IF (kflag<0) THEN
        namest(ncall) = 'FMMPYE'
        CALL fmwarn
      END IF

      IF (ma2*mb2<0) me(2) = -me(2)
      IF (ma2*mc2<0) mf(2) = -mf(2)
      IF (ma2*md2<0) mg(2) = -mg(2)

      IF (kaccsw==1) THEN
        md2b = nint((ndig-1)*alogm2+log(real(abs(me(2))+1))/0.69315)
        me(0) = min(macca,maccb,md2b)
        md2b = nint((ndig-1)*alogm2+log(real(abs(mf(2))+1))/0.69315)
        mf(0) = min(macca,maccc,md2b)
        md2b = nint((ndig-1)*alogm2+log(real(abs(mg(2))+1))/0.69315)
        mg(0) = min(macca,maccd,md2b)
      ELSE
        me(0) = min(macca,maccb)
        mf(0) = min(macca,maccc)
        mg(0) = min(macca,maccd)
      END IF

120   IF (ntrace/=0) THEN
        CALL fmntr(1,me,me,1)
        IF (abs(ntrace)>=1 .AND. ncall<=lvltrc) THEN
          IF (ntrace<0) THEN
            CALL fmntrj(mf,ndig)
            CALL fmntrj(mg,ndig)
          ELSE
            CALL fmprnt(mf)
            CALL fmprnt(mg)
          END IF
        END IF
      END IF
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmmpye
    SUBROUTINE fmmpyi(ma,ival,mb)

!  MB = MA * IVAL

!  Multiply FM number MA by one word integer IVAL.

!  This routine is faster than FMMPY when IVAL*MBASE is a
!  one word integer.

      IMPLICIT NONE

!             Scratch array usage during FMMPYI:   M01

! .. Intrinsic Functions ..
      INTRINSIC abs, dble, int, log, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, macca, mcarry, md2b, mkt, mlr, mval
      INTEGER :: j, ka, kb, kc, kshift, n1, nguard, nmval, nv2
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmeq, fmim, fmmove, fmmpy2, fmntr, fmntri, fmrnd, &
        fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      macca = ma(0)
      ncall = ncall + 1
      IF (ntrace/=0) THEN
        namest(ncall) = 'FMMPYI'
        CALL fmntr(2,ma,ma,1)
        CALL fmntri(2,ival,0)
      END IF
      kflag = 0
      n1 = ndig + 1

!             Check for special cases.

      IF (ma(2)==0) THEN
        CALL fmeq(ma,mb)
        IF (ntrace/=0) THEN
          CALL fmntr(1,mb,mb,1)
        END IF
        ncall = ncall - 1
        RETURN
      END IF

      IF (abs(ma(1))<mexpov .AND. abs(ival)>1) GO TO 20

      IF (ma(1)==munkno) THEN
        CALL fmim(0,mb)
        mb(1) = munkno
        mb(2) = 1
        mb(0) = nint(ndig*alogm2)
        kflag = -4
        IF (ntrace/=0) THEN
          CALL fmntr(1,mb,mb,1)
        END IF
        ncall = ncall - 1
        RETURN
      END IF

      IF (ival==0) THEN
        CALL fmim(0,mb)
        IF (ntrace/=0) THEN
          CALL fmntr(1,mb,mb,1)
        END IF
        ncall = ncall - 1
        RETURN
      END IF

      IF (abs(ival)==1) THEN
        DO 10 j = 0, n1
          mb(j) = ma(j)
10      CONTINUE
        IF (ma(1)==mexpov) kflag = -5
        IF (ma(1)==mexpun) kflag = -6
        mb(2) = ma(2)*ival
        IF (ntrace/=0) THEN
          CALL fmntr(1,mb,mb,1)
        END IF
        ncall = ncall - 1
        RETURN
      END IF

      IF (ma(1)==mexpov) THEN
        ma2 = ma(2)
        CALL fmim(0,mb)
        kflag = -5
        mb(1) = mexpov
        mb(2) = 1
        mb(0) = nint(ndig*alogm2)
        IF ((ma2<0 .AND. ival>0) .OR. (ma2>0 .AND. ival<0)) mb(2) = -1
        IF (ntrace/=0) THEN
          CALL fmntr(1,mb,mb,1)
        END IF
        ncall = ncall - 1
        RETURN
      END IF

      IF (ma(1)==mexpun) THEN
        CALL fmim(0,mb)
        mb(1) = munkno
        mb(2) = 1
        mb(0) = nint(ndig*alogm2)
        namest(ncall) = 'FMMPYI'
        kflag = -4
        CALL fmwarn
        IF (ntrace/=0) THEN
          CALL fmntr(1,mb,mb,1)
        END IF
        ncall = ncall - 1
        RETURN
      END IF

!             Work with positive numbers.

20    ma2 = ma(2)
      ma(2) = abs(ma(2))
      mval = abs(ival)
      nmval = int(mval)
      nv2 = nmval - 1

!             To leave room for the normalization, shift the product
!             to the right KSHIFT places in MWA.

      kshift = int((log(dble(ma(2)+1)*dble(mval)))/dlogmb)

!             If IVAL is too big use FMMPY.

      IF (kshift>ndig .OR. mval>maxint/mbase .OR. nmval/=abs(ival) .OR. &
          nv2/=abs(ival)-1) THEN
        CALL fmim(ival,m01)
        ma(2) = ma2
        CALL fmmpy2(ma,m01,mb)
        IF (ntrace/=0) THEN
          CALL fmntr(1,mb,mb,1)
        END IF
        ncall = ncall - 1
        RETURN
      END IF

      mwa(1) = ma(1) + kshift
      ka = 2 + kshift
      kb = n1 + kshift
      kc = ndig + 5
      DO 30 j = kb, kc
        mwa(j) = 0
30    CONTINUE

      mcarry = 0

!             This is the main multiplication loop.

      DO 40 j = kb, ka, -1
        mkt = ma(j-kshift)*mval + mcarry
        mcarry = int(mkt/mbase)
        mwa(j) = mkt - mcarry*mbase
40    CONTINUE

!             Resolve the final carry.

      DO 50 j = ka - 1, 2, -1
        mkt = int(mcarry/mbase)
        mwa(j) = mcarry - mkt*mbase
        mcarry = mkt
50    CONTINUE

!             Now the first significant digit in the product is in
!             MWA(2) or MWA(3).  Round the result and move it to MB.

      ma(2) = ma2
      IF (mwa(2)==0) THEN
        mlr = 2*mwa(ndig+3) + 1
        IF (mlr>=mbase) THEN
          IF (mlr-1>mbase .AND. mwa(n1+1)<mbase-1) THEN
            IF (kround/=0 .OR. ncall>1) THEN
              mwa(n1+1) = mwa(n1+1) + 1
              mwa(n1+2) = 0
            END IF
          ELSE
            nguard = kshift - 1
            CALL fmrnd(mwa,ndig,nguard,1)
          END IF
        END IF
      ELSE
        mlr = 2*mwa(ndig+2) + 1
        IF (mlr>=mbase) THEN
          IF (mlr-1>mbase .AND. mwa(n1)<mbase-1) THEN
            IF (kround/=0 .OR. ncall>1) THEN
              mwa(n1) = mwa(n1) + 1
              mwa(n1+1) = 0
            END IF
          ELSE
            CALL fmrnd(mwa,ndig,kshift,0)
          END IF
        END IF
      END IF
      CALL fmmove(mwa,mb)

      IF (kflag<0) THEN
        namest(ncall) = 'FMMPYI'
        CALL fmwarn
      END IF

!             Put the sign on the result.

      IF ((ival>0 .AND. ma2<0) .OR. (ival<0 .AND. ma2>0)) mb(2) = -mb(2)

      IF (kaccsw==1) THEN
        md2b = nint((ndig-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
        mb(0) = min(macca,md2b)
      ELSE
        mb(0) = macca
      END IF

      IF (ntrace/=0) THEN
        CALL fmntr(1,mb,mb,1)
      END IF
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmmpyi
    SUBROUTINE fmmset(maxint,ml,mld2,mlm1)

!  Internal routine to keep some compilers from doing a loop at
!  the highest precision available and then rounding to the
!  declared precision.  For example, it is used in FMSET while
!  trying to find the threshold beyond which integers cannot
!  be represented exactly using (M) precision.

! .. Intrinsic Functions ..
      INTRINSIC dint
! ..
! .. Scalar Arguments ..
      REAL (KIND(0.0D0)) :: maxint, ml, mld2, mlm1
! ..
      ml = 2*maxint + 1
      mld2 = dint(ml/2)
      mlm1 = ml - 1

      RETURN
    END SUBROUTINE fmmset
    SUBROUTINE fmnint(ma,mb)

!  MB = NINT(MA)  --  MB is returned as the nearest integer to MA.

      IMPLICIT NONE

!             Scratch array usage during FMNINT:   M01

! .. Intrinsic Functions ..
      INTRINSIC abs, int, max
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, mxsave
      INTEGER :: k, kasave, kovun, kreslt, kwrnsv, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmdivi, fmentr, fmeq2, fmexit, fmi2m, fmint, &
        fmntr, fmrslt, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (abs(ma(1))>mexpab) THEN
        CALL fmentr('FMNINT',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        namest(ncall) = 'FMNINT'
        IF (ntrace/=0) CALL fmntr(2,ma,ma,1)
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
        ndsave = ndig
        IF (ncall==1) THEN
          k = max(ngrd52-1,2)
          ndig = max(ndig+k,2)
          IF (ndig>ndg2mx) THEN
            kflag = -9
            CALL fmwarn
            ndig = ndsave
            kreslt = 12
            CALL fmrslt(ma,ma,mb,kreslt)
            IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
            ncall = ncall - 1
            RETURN
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 0
        mxsave = mxexp
        mxexp = mxexp2
      END IF

      kwrnsv = kwarn
      kwarn = 0
      CALL fmeq2(ma,mb,ndsave,ndig,0)
      IF (ndsave>int(ma(1))) THEN
        ma2 = ma(2)
        mb(2) = abs(mb(2))
        CALL fmi2m(1,m01)
        CALL fmdivi(m01,2,m01)
        CALL fmadd(mb,m01,mb)
        CALL fmint(mb,mb)
        IF (ma2<0) mb(2) = -mb(2)
      END IF
      kwarn = kwrnsv

!             Round the result and return.

      CALL fmexit(mb,mb,ndsave,mxsave,kasave,kovun)
      RETURN
    END SUBROUTINE fmnint
    SUBROUTINE fmntr(ntr,ma,mb,narg)

!  Print FM numbers in base 10 format using FMOUT for conversion.
!  This is used for trace output from the FM routines.

!  NTR =  1 if a result of an FM call is to be printed.
!      =  2 to print input argument(s) to an FM call.

!  MA  -  the FM number to be printed.

!  MB  -  an optional second FM number to be printed.

!  NARG - the number of arguments.  NARG = 1 if only MA is to be
!         printed, and NARG = 2 if both MA and MB are to be printed.

!  NTRACE and LVLTRC (in COMMON /FMUSER/) control trace printout.

!  NTRACE = 0        No printout except warnings and errors.

!  NTRACE = 1        The result of each call to one of the routines
!                    is printed in base 10, using FMOUT.

!  NTRACE = -1       The result of each call to one of the routines
!                    is printed in internal base MBASE format.

!  NTRACE = 2        The input arguments and result of each call to one
!                    of the routines is printed in base 10, using FMOUT.

!  NTRACE = -2       The input arguments and result of each call to one
!                    of the routines is printed in base MBASE format.

!  LVLTRC defines the call level to which the trace is done.  LVLTRC = 1
!         means only FM routines called directly by the user are traced,
!         LVLTRC = K prints traces for FM routines with call levels up
!         to and including level K.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, int
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: narg, ntr
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      CHARACTER (6) :: name
! ..
! .. External Subroutines ..
      EXTERNAL fmntrj, fmprnt
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (ntrace==0) RETURN
      IF (ncall>lvltrc) RETURN
      IF (ntr==2 .AND. abs(ntrace)==1) RETURN

      IF (ntr==2) THEN
        name = namest(ncall)
        WRITE (kw,90000) name
      ELSE
        name = namest(ncall)
        IF (kflag==0) THEN
          WRITE (kw,90010) name, ncall, int(mbase), ndig
        ELSE
          WRITE (kw,90020) name, ncall, int(mbase), ndig, kflag
        END IF
      END IF

!             Check for base MBASE internal format trace.

      IF (ntrace<0) THEN
        CALL fmntrj(ma,ndig)
        IF (narg==2) CALL fmntrj(mb,ndig)
      END IF

!             Check for base 10 trace using FMOUT.

      IF (ntrace>0) THEN
        CALL fmprnt(ma)

        IF (narg==2) THEN
          CALL fmprnt(mb)
        END IF
      END IF

      RETURN
90000 FORMAT (' Input to ',A6)
90010 FORMAT (' ',A6,15X,'Call level =',I2,5X,'MBASE =',I10,5X,'NDIG =',I6)
90020 FORMAT (' ',A6,6X,'Call level =',I2,4X,'MBASE =',I10,4X,'NDIG =',I6,4X, &
        'KFLAG =',I3)
    END SUBROUTINE fmntr
    SUBROUTINE fmntri(ntr,n,knam)

!  Internal routine for trace output of integer variables.

!  NTR = 1 for output values
!        2 for input values

!  N     Integer to be printed.

!  KNAM  is positive if the routine name is to be printed.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, int
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: knam, n, ntr
! ..
! .. Local Scalars ..
      CHARACTER (6) :: name
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (ntrace==0) RETURN
      IF (ncall>lvltrc) RETURN
      IF (ntr==2 .AND. abs(ntrace)==1) RETURN

      IF (ntr==2 .AND. knam>0) THEN
        name = namest(ncall)
        WRITE (kw,90000) name
      END IF
      IF (ntr==1 .AND. knam>0) THEN
        name = namest(ncall)
        IF (kflag==0) THEN
          WRITE (kw,90010) name, ncall, int(mbase), ndig
        ELSE
          WRITE (kw,90020) name, ncall, int(mbase), ndig, kflag
        END IF
      END IF

      WRITE (kw,90030) n

      RETURN
90000 FORMAT (' Input to ',A6)
90010 FORMAT (' ',A6,15X,'Call level =',I2,5X,'MBASE =',I10,5X,'NDIG =',I6)
90020 FORMAT (' ',A6,6X,'Call level =',I2,4X,'MBASE =',I10,4X,'NDIG =',I6,4X, &
        'KFLAG =',I3)
90030 FORMAT (1X,I18)
    END SUBROUTINE fmntri
    SUBROUTINE fmntrj(ma,nd)

!  Print trace output in internal base MBASE format.  The number to
!  be printed is in MA.

!  ND is the number of base MBASE digits to be printed.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC dble, int, log10
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: nd
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: j, l, n, n1
      CHARACTER (50) :: form
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      n1 = nd + 1

      l = int(log10(dble(mbase-1))) + 2
      n = (kswide-23)/l
      IF (n>10) n = 5*(n/5)
      IF (nd<=n) THEN
        WRITE (form,90000) l + 2, n - 1, l
      ELSE
        WRITE (form,90010) l + 2, n - 1, l, n, l
      END IF
      WRITE (kw,form) (int(ma(j)),j=1,n1)

      RETURN
90000 FORMAT (' (1X,I19,I',I2,',',I3,'I',I2,') ')
90010 FORMAT (' (1X,I19,I',I2,',',I3,'I',I2,'/(22X,',I3,'I',I2,')) ')
    END SUBROUTINE fmntrj
    SUBROUTINE fmntrr(ntr,x,knam)

!  Internal routine for trace output of real variables.

!  NTR - 1 for output values
!        2 for input values

!  X   - Double precision value to be printed if NX.EQ.1

!  KNAM - Positive if the routine name is to be printed.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, int
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL (KIND(0.0D0)) :: x
      INTEGER :: knam, ntr
! ..
! .. Local Scalars ..
      CHARACTER (6) :: name
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (ntrace==0) RETURN
      IF (ncall>lvltrc) RETURN
      IF (ntr==2 .AND. abs(ntrace)==1) RETURN

      IF (ntr==2 .AND. knam>0) THEN
        name = namest(ncall)
        WRITE (kw,90000) name
      END IF
      IF (ntr==1 .AND. knam>0) THEN
        name = namest(ncall)
        IF (kflag==0) THEN
          WRITE (kw,90010) name, ncall, int(mbase), ndig
        ELSE
          WRITE (kw,90020) name, ncall, int(mbase), ndig, kflag
        END IF
      END IF

      WRITE (kw,90030) x

      RETURN
90000 FORMAT (' Input to ',A6)
90010 FORMAT (' ',A6,15X,'Call level =',I2,5X,'MBASE =',I10,5X,'NDIG =',I6)
90020 FORMAT (' ',A6,6X,'Call level =',I2,4X,'MBASE =',I10,4X,'NDIG =',I6,4X, &
        'KFLAG =',I3)
90030 FORMAT (1X,D30.20)
    END SUBROUTINE fmntrr
    SUBROUTINE fmout(ma,line,lb)

!  Convert a floating multiple precision number to a character array
!  for output.

!  MA   is an FM number to be converted to an A1 character
!       array in base 10 format
!  LINE is the CHARACTER*1 array in which the result is returned.
!  LB   is the length of LINE.

! JFORM1 and JFORM2 (in COMMON) determine the format of LINE.

! JFORM1 = 0  normal setting  ( .314159M+6 )
!        = 1  1PE format      ( 3.14159M+5 )
!        = 2  F   format      ( 314159.000 )

! JFORM2 = number of significant digits to display (if JFORM1 = 0, 1)
!        = number of digits after the decimal point (if JFORM1 = 2)

!          If JFORM2.EQ.0 and JFORM1.NE.2 then a default number of
!          digits is chosen.  The default is roughly the full precision
!          of MA.

!          If JFORM2.EQ.0 and JFORM1.EQ.2 then the number is returned in
!          integer format with no decimal point.  Rounding is done as
!          with other settings, so the value displayed is the nearest
!          integer to MA.

!  If JFORM1.EQ.2 and MA is too large or too small to display in the
!  requested format, it is converted using JFORM1=0, JFORM2=0.

!  LINE should be dimensioned at least LOG10(MBASE)*NDIG + 15 on a
!  32-bit machine to allow for up to 10 digit exponents.  Replace
!  15 by 20 if 48-bit integers are used, 25 for 64-bit integers, ....

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dble, dint, int, log10, max, min, mod, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: lb
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
      CHARACTER (1) :: line(lb)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: m2, mbsave, mexp, mexp10, mkt, mndgms, ms1, ms2, msd2, mt10, &
        mxsave
      REAL :: x
      INTEGER :: j, jdpt, jf1sav, jf2sav, k, k1, k2, ka, kasave, kb, kc, &
        kdigit, kexp, kexpsh, kms2sd, kmt, kpt, krsave, l, nd, nde, nde2, &
        ndigms, nds2, ndsave, npower, nsd1, nsd2, nval, nword, nword1, nword2
      CHARACTER (1) :: kchar
! ..
! .. Local Arrays ..
      REAL (KIND(0.0D0)) :: md(0:lunpck), ms(0:lunpck), mt(0:lunpck)
      CHARACTER (1) :: nexpov(12), nexpun(12), numb(10), nunkno(12)
! ..
! .. External Subroutines ..
      EXTERNAL fmadd2, fmcons, fmdiv2, fmeq, fmeq2, fmim, fmmpy2, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
! .. Data Statements ..
      DATA numb/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/
      DATA nunkno/' ', ' ', ' ', 'U', 'N', 'K', 'N', 'O', 'W', 'N', ' ', ' '/
      DATA nexpov/' ', ' ', ' ', 'O', 'V', 'E', 'R', 'F', 'L', 'O', 'W', ' '/
      DATA nexpun/' ', ' ', ' ', 'U', 'N', 'D', 'E', 'R', 'F', 'L', 'O', 'W'/
! ..
!             To avoid recursion, FMOUT calls only internal arithmetic
!             routines (FMADD2, FMMPY2, ...), so no trace printout is
!             done during a call to FMOUT.

      kflag = 0
      ncall = ncall + 1
      namest(ncall) = 'FMOUT '

!             Raise the call stack again, since the internal
!             routines don't.

      ncall = ncall + 1
      namest(ncall) = 'FMOUT '
      DO 10 j = 1, lb
        line(j) = ' '
10    CONTINUE

!             Check for special cases.

      IF (ma(1)==munkno) THEN
        DO 20 j = 1, 12
          line(j) = nunkno(j)
20      CONTINUE
        ncall = ncall - 2
        RETURN
      END IF
      IF (ma(1)==mexpov) THEN
        DO 30 j = 1, 12
          line(j) = nexpov(j)
30      CONTINUE
        line(2) = '+'
        IF (ma(2)<0) line(2) = '-'
        ncall = ncall - 2
        RETURN
      END IF
      IF (ma(1)==mexpun) THEN
        DO 40 j = 1, 12
          line(j) = nexpun(j)
40      CONTINUE
        line(2) = '+'
        IF (ma(2)<0) line(2) = '-'
        ncall = ncall - 2
        RETURN
      END IF
      IF (ma(2)==0 .AND. jform1==2 .AND. jform2==0) THEN
        line(2) = '0'
        ncall = ncall - 2
        RETURN
      END IF

      kasave = kaccsw
      kaccsw = 0
      krsave = kround
      kround = 1
      jf1sav = jform1
      jf2sav = jform2
      mbsave = mbase
      ndsave = ndig
      mxsave = mxexp

!             ND is the number of base 10 digits required.

50    nd = jform2
      IF (jform1==2 .AND. ma(1)>0) nd = jform2 + int(real(ma(1))*log10(real( &
        mbase))) + 1
      IF (nd<=1) THEN
        k = int(real(ndig)*log10(real(mbase)))
        nd = max(k,jform2)
      END IF
      IF (jform2<=0 .AND. jform1<=1) nd = int(1.1+real(ndig-1)*log10(real( &
        mbase)))
      IF (nd<2) nd = 2

      IF (lb<nd+6) THEN
        IF (jform1==2) THEN
          jform1 = 0
          jform2 = 0
          GO TO 50
        END IF
        GO TO 270
      END IF

!             Convert to the base that is the largest power of 10
!             less than MXBASE and build the output number.

      npower = int(log10(real(mxbase)/4))
      mxexp = mxexp2
      mbase = 10**npower
      IF (mblogs/=mbase) CALL fmcons
      ndig = nd/npower + 3
      IF (ndig<2) ndig = 2
      IF (ndig>ndg2mx) THEN
        kflag = -9
        ncall = ncall - 1
        CALL fmwarn
        ncall = ncall + 1
        GO TO 270
      END IF

      IF (ma(2)==0) THEN
        CALL fmim(0,ms)
        GO TO 110
      END IF

!             Check to see if MA is already in a base that is a
!             power of ten.  If so, the conversion can be skipped.

      k = npower
      DO 60 j = 1, k
        mbase = 10**j
        IF (mbase==mbsave) THEN
          IF (mblogs/=mbase) CALL fmcons
          npower = j
          ndig = nd/npower + 2
          IF (ndig<2) ndig = 2
          IF (ndig>ndg2mx) THEN
            kflag = -9
            ncall = ncall - 1
            CALL fmwarn
            ncall = ncall + 1
            GO TO 270
          END IF
          CALL fmeq2(ma,ms,ndsave,ndig,0)
          ms(2) = abs(ms(2))
          GO TO 110
        END IF
60    CONTINUE

      IF (mblogs/=mbase) CALL fmcons
      CALL fmim(int(mbsave),md)
      nds2 = ndsave + 1
      CALL fmim(1,mt)
      kmt = 1

!             Convert the fraction part of MA to the new base.

      kpt = nds2 + 1
      DO 70 j = 3, nds2
        kpt = kpt - 1
        IF (ma(kpt)/=0) GO TO 80
70    CONTINUE

80    kexpsh = kpt - 1
      kdigit = int(abs(ma(2)))
      CALL fmim(kdigit,ms)
      ndigms = ndig

      DO 90 j = 3, kpt
        kdigit = int(ma(j))
        IF (mbsave==2) THEN
          ndig = min(ndigms,max(2,int(ms(1))+1))
          CALL fmadd2(ms,ms,ms)
        ELSE
          ndig = min(ndigms,max(2,int(ms(1)+md(1))))
          CALL fmmpy2(ms,md,ms)
        END IF

        IF (kdigit>0) THEN
          IF (kmt/=kdigit) THEN
            ndig = min(ndigms,max(2,int(md(1))))
            CALL fmim(kdigit,mt)
            kmt = kdigit
          END IF
          ndig = min(ndigms,max(2,int(max(ms(1),mt(1)))+1))
          CALL fmadd2(ms,mt,ms)
        END IF
90    CONTINUE

!             Convert the exponent.

      ndig = ndigms
      CALL fmim(1,mt)
      k = abs(int(ma(1))-kexpsh)
      IF (mod(k,2)==1) THEN
        CALL fmeq(md,mt)
      ELSE
        CALL fmim(1,mt)
      END IF

100   k = k/2
      m2 = 2
      mndgms = ndigms
      ndig = int(min(mndgms,max(m2,md(1)*m2)))
      IF (k>0) CALL fmmpy2(md,md,md)
      IF (mod(k,2)==1) THEN
        ndig = int(min(mndgms,max(m2,mt(1)+md(1))))
        CALL fmmpy2(mt,md,mt)
      END IF
      IF (k>1) GO TO 100

      ndig = ndigms
      IF (ma(1)-kexpsh<0) THEN
        CALL fmdiv2(ms,mt,ms)
      ELSE
        CALL fmmpy2(ms,mt,ms)
      END IF

!             Now MS is the value of MA converted to a
!             power of ten base.

!             Convert it to a character string base 10 for output.

!             MEXP10 is the base 10 exponent.
!             KMS2SD is the number of base 10 significant digits
!                    in MS(2).

110   ms1 = ms(1)
120   mexp10 = npower*ms(1)
      kms2sd = npower
      k = int(mbase)
      DO 130 j = 1, npower
        k = k/10
        IF (ms(2)<k .AND. ms(2)/=0) THEN
          mexp10 = mexp10 - 1
          kms2sd = kms2sd - 1
        END IF
130   CONTINUE

!             For printing using JFORM1 = 1, reduce the exponent to
!             account for the fact that the decimal point and first
!             significant digit will later be swapped.

      IF (jform1==1 .AND. ms(2)/=0) mexp10 = mexp10 - 1

!             Find the position in the unpacked number for rounding.
!             NWORD is the word in which rounding is done, or zero if
!                   no rounding is necessary.
!                   NWORD is set to -1 if JFORM1 is 2 (F format) but no
!                   significant digits would be printed.  This case
!                   defaults to JFORM1 = 0.
!             NVAL gives the position within that word where rounding
!                  occurs.
!             NSD1 is the maximum number of base 10 S.D.'s in NWORD
!                  digits of base 10**NPOWER.
!             NSD2 is the number of base 10 S.D.'s needed to get ND
!                  base 10 digits after the decimal.

      nsd2 = nd
      IF (jform1==2) THEN
        msd2 = jform2 + mexp10
        IF (msd2>nd) THEN
          nsd2 = nd
        ELSE
          nsd2 = int(msd2)
        END IF
        nword = (nsd2-kms2sd-1+npower)/npower + 2
        IF (nword<2) nword = -1
        IF (nword>ndig) nword = 0
        IF (nword>=2 .AND. nsd2<=0) nword = -1
      ELSE
        nword = (nd-kms2sd-1+npower)/npower + 2
      END IF
      nsd1 = kms2sd + npower*(nword-2)
      IF (nword<2) THEN
        nval = 0
      ELSE
        nval = 10**(nsd1-nsd2)
      END IF

!             Now do the base 10 rounding.

      IF (nword>=2) THEN
        x = 0.0
        IF (nval>1) x = mod(int(ms(nword)),nval)
        IF (nword<ndig+1) THEN
          x = real(dble(x)+dble(ms(nword+1))/dble(mbase))
        END IF
        x = x/nval
        IF (x<0.5) GO TO 140
        ms2 = ms(2)
        ms(nword) = int(ms(nword)/nval)*nval
        ms(nword+1) = 0
        ms(nword+2) = 0
        ms(nword) = ms(nword) + nval
        IF (ms(nword)>=mbase) THEN
          nword1 = nword - 1
          nword2 = nword - 2
          IF (nword>2) THEN
            CALL fmeq2(ms,ms,nword1,nword2,1)
          ELSE
            ms(1) = ms(1) + 1
            ms(2) = int(ms(2)/mbase)
            ms(3) = 0
          END IF
        END IF
        IF (ms(1)/=ms1 .OR. ms(2)/=ms2) GO TO 120
      END IF

!             Build the base 10 character string.

140   IF (ma(2)<0) line(1) = '-'
      line(2) = '.'
      k = 10**kms2sd
      l = 2
      IF (nword==-1) nsd2 = nd
      DO 150 j = 1, nsd2
        k = k/10
        IF (k==0) THEN
          k = int(mbase)/10
          l = l + 1
        END IF
        kdigit = int(ms(l))/k
        ms(l) = mod(int(ms(l)),k)
        line(j+2) = numb(kdigit+1)
150   CONTINUE

      ka = nsd2 + 3
      kb = nd + 2
      IF (kb>=ka) THEN
        DO 160 j = ka, kb
          line(j) = numb(1)
160     CONTINUE
      END IF

      line(nd+3) = cmchar
      line(nd+4) = '+'
      IF (mexp10<0) line(nd+4) = '-'
      IF (ma(2)==0) line(nd+4) = ' '

!             Build the digits of the base 10 exponent backwards,
!             then reverse them.

      nde = 1
      mexp = abs(mexp10)
      mt10 = 10
      DO 180 j = 1, lb
        mkt = dint(mexp/mt10)
        kdigit = int(mexp-mkt*mt10)
        line(nd+4+j) = numb(kdigit+1)
        mexp = mkt
        IF (mexp==0) GO TO 190

        IF (nd+5+j>lb) THEN
          DO 170 k = 1, lb
            line(k) = '*'
170       CONTINUE
          GO TO 210
        END IF

        nde = nde + 1
180   CONTINUE

190   nde2 = nde/2
      IF (nde2<1) GO TO 210
      k1 = nd + 4
      k2 = nd + 5 + nde
      DO 200 j = 1, nde2
        k1 = k1 + 1
        k2 = k2 - 1
        kchar = line(k1)
        line(k1) = line(k2)
        line(k2) = kchar
200   CONTINUE

!             If JFORM1 is 1 put the first digit left of the decimal.

210   IF (jform1==1) THEN
        kchar = line(2)
        line(2) = line(3)
        line(3) = kchar
      END IF

!             If JFORM1 is 2 put the number into fixed format.

      IF (jform1==2 .AND. jform2>=0) THEN
        IF (mexp10<=-jform2 .OR. mexp10+2>lb) THEN
          jform1 = 0
          jform2 = 0
          mbase = mbsave
          IF (mblogs/=mbase) CALL fmcons
          ndig = ndsave
          mxexp = mxsave
          DO 220 j = 1, lb
            line(j) = ' '
220       CONTINUE
          GO TO 50
        END IF
        ka = nd + 3
        DO 230 j = ka, lb
          line(j) = numb(1)
230     CONTINUE

        kexp = int(mexp10)
        IF (mexp10>0) THEN
          DO 240 j = 1, kexp
            line(j+1) = line(j+2)
240       CONTINUE
          line(kexp+2) = '.'
        END IF

        IF (mexp10<0) THEN
          kexp = -int(mexp10)
          ka = 3 + kexp
          kb = lb + 1
          kc = kb - kexp
          DO 250 j = ka, lb
            kb = kb - 1
            kc = kc - 1
            line(kb) = line(kc)
            line(kc) = numb(1)
250       CONTINUE
        END IF

        jdpt = 0
        DO 260 j = 1, lb
          IF (line(j)=='.') jdpt = j
          IF (jdpt>0 .AND. j>jdpt+jform2) line(j) = ' '
260     CONTINUE
        IF (jform2==0 .AND. jdpt>0) line(kexp+2) = ' '

      END IF

!             Restore values and return

      GO TO 290

!             LINE is not big enough to hold the number
!             of digits specified.

270   kflag = -8
      DO 280 j = 1, lb
        line(j) = '*'
280   CONTINUE
      ncall = ncall - 1
      CALL fmwarn
      ncall = ncall + 1

290   mbase = mbsave
      IF (mblogs/=mbase) CALL fmcons
      ndig = ndsave
      mxexp = mxsave
      ncall = ncall - 2
      kaccsw = kasave
      kround = krsave
      jform1 = jf1sav
      jform2 = jf2sav
      RETURN
    END SUBROUTINE fmout
    SUBROUTINE fmpack(ma,mp)

!  MA is packed two base NDIG digits per word and returned in MP.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, mod
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mp(0:lpack)
! ..
! .. Local Scalars ..
      INTEGER :: j, kp
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      kp = 2
      mp(0) = ma(0)
      mp(1) = ma(1)
      mp(2) = abs(ma(2))*mbase + ma(3)
      IF (ma(2)<0) mp(2) = -mp(2)
      IF (ndig>=4) THEN
        DO 10 j = 4, ndig, 2
          kp = kp + 1
          mp(kp) = ma(j)*mbase + ma(j+1)
10      CONTINUE
      END IF
      IF (mod(ndig,2)==1) mp(kp+1) = ma(ndig+1)*mbase
      RETURN
    END SUBROUTINE fmpack
    SUBROUTINE fmpi(ma)

!  MA = pi

      IMPLICIT NONE

!             Scratch array usage during FMPI:   M01 - M04

! .. Intrinsic Functions ..
      INTRINSIC abs, int, max, min, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: j, k, kasave, ndmb, ndsave, ndsv
      CHARACTER (155) :: string
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmeq2, fmntr, fmpi2, fmst2m, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      kflag = 0
      ncall = ncall + 1
      namest(ncall) = 'FMPI  '
      IF (abs(ntrace)>=2 .AND. ncall<=lvltrc) THEN
        WRITE (kw,90000)
      END IF
      kasave = kaccsw
      kaccsw = 0

!             Increase the working precision.

      ndsave = ndig
      IF (ncall==1) THEN
        k = ngrd52
        ndig = max(ndig+k,2)
        IF (ndig>ndg2mx) THEN
          kflag = -9
          CALL fmwarn
          ma(1) = munkno
          ma(2) = 1
          ma(0) = nint(ndig*alogm2)
          DO 10 j = 2, ndsave
            ma(j+1) = 0
10        CONTINUE
          GO TO 20
        END IF
      END IF

!             Check to see if pi has previously been computed
!             in base MBASE with sufficient precision.

      IF (mbspi==mbase .AND. ndigpi>=ndig) THEN
        IF (namest(ncall-1)/='NOEQ  ') THEN
          kaccsw = kasave
          CALL fmeq2(mpisav,ma,ndigpi,ndsave,0)
        END IF
      ELSE
        ndmb = int(150.0*2.302585/alogmb)
        IF (ndmb>=ndig) THEN
          ndsv = ndig
          ndig = min(ndmb,ndg2mx)
          string = '3.141592653589793238462643383279502884197169' // &
            '39937510582097494459230781640628620899862803482534211' // &
            '7067982148086513282306647093844609550582231725359408128'
          CALL fmst2m(string,mpisav)
          mpisav(0) = nint(ndig*alogm2)
          mbspi = mbase
          ndigpi = ndig
          IF (abs(mpisav(1))>10) ndigpi = 0
        ELSE
          ndsv = ndig
          ndig = min(ndig+2,ndg2mx)
          CALL fmpi2(mpisav)
          mpisav(0) = nint(ndig*alogm2)
          mbspi = mbase
          ndigpi = ndig
          IF (abs(mpisav(1))>10) ndigpi = 0
        END IF
        IF (namest(ncall-1)/='NOEQ  ') THEN
          kaccsw = kasave
          CALL fmeq2(mpisav,ma,ndig,ndsave,0)
        END IF
        ndig = ndsv
      END IF

20    ndig = ndsave
      kaccsw = kasave
      IF (ntrace/=0) CALL fmntr(1,ma,ma,1)
      ncall = ncall - 1
      RETURN
90000 FORMAT (' Input to FMPI')
    END SUBROUTINE fmpi
    SUBROUTINE fmpi2(mpi)

!  Internal routine to compute pi.
!  The formula used is due to S. Ramanujan:
!                                                (4n)!(1103+26390n)
!  1/pi = (sqrt(8)/9801) * sum(n=0 to infinity) --------------------
!                                               ((n!)**4)(396**(4n))
!  The result is returned in MPI.

      IMPLICIT NONE

!             Scratch array usage during FMPI2:   M01 - M04

! .. Intrinsic Functions ..
      INTRINSIC int, max, nint, sqrt
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: mpi(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: mx
      REAL (KIND(0.0D0)) :: x
      INTEGER :: j, k, kst, large, n, ndigrd, ndsave
! ..
! .. Local Arrays ..
      INTEGER :: nstack(19)
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmdig, fmdiv, fmdivi, fmdpm, fmi2m, fmmpy, &
        fmmpyi
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      ndsave = ndig
      n = -1
      CALL fmi2m(1103,mpi)
      CALL fmi2m(1,m02)
      CALL fmi2m(26390,m03)
      CALL fmi2m(1103,m04)
      mx = mxbase**2/mbase
      IF (mx>mxexp2) mx = mxexp2

10    n = n + 1
      large = int(mx)/(4*n+3)
      j = 4*n + 1
      IF (j>large) THEN
        CALL fmmpyi(m02,j,m02)
        j = j + 1
        CALL fmmpyi(m02,j,m02)
        j = j + 1
        CALL fmmpyi(m02,j,m02)
      ELSE IF (j*(j+1)>large) THEN
        k = j*(j+1)
        CALL fmmpyi(m02,k,m02)
        j = j + 2
        CALL fmmpyi(m02,j,m02)
      ELSE
        k = j*(j+1)*(j+2)
        CALL fmmpyi(m02,k,m02)
      END IF

      j = n + 1
      large = int(mxbase)/j
      IF (j>large) THEN
        CALL fmdivi(m02,j,m02)
        CALL fmdivi(m02,j,m02)
        CALL fmdivi(m02,j,m02)
      ELSE IF (j*j>large) THEN
        k = j*j
        CALL fmdivi(m02,k,m02)
        CALL fmdivi(m02,j,m02)
      ELSE
        k = j*j*j
        CALL fmdivi(m02,k,m02)
      END IF

!             Break 4/396**4 into 1/(2178*2178*1296).

      j = 2178
      large = int(mxbase)/j
      IF (j>large) THEN
        CALL fmdivi(m02,j,m02)
        CALL fmdivi(m02,j,m02)
        CALL fmdivi(m02,1296,m02)
      ELSE
        k = j*j
        CALL fmdivi(m02,k,m02)
        CALL fmdivi(m02,1296,m02)
      END IF

      ndigrd = ndig
      ndig = ndsave
      CALL fmadd(m03,m04,m04)
      ndig = ndigrd
      CALL fmmpy(m02,m04,m01)

      ndig = ndsave
      CALL fmadd(mpi,m01,mpi)
      ndig = max(2,ndsave-int(mpi(1)-m01(1)))
      IF (kflag/=1) GO TO 10
      ndig = ndsave

      CALL fmi2m(8,m02)
      x = 8
      x = sqrt(x)
      CALL fmdpm(x,m04)
      CALL fmdig(nstack,kst)
      DO 20 j = 1, kst
        ndig = nstack(j)
        CALL fmdiv(m02,m04,m01)
        CALL fmadd(m04,m01,m04)
        CALL fmdivi(m04,2,m04)
20    CONTINUE
      m04(0) = nint(ndig*alogm2)
      CALL fmi2m(9801,m03)
      CALL fmmpy(mpi,m04,mpi)
      CALL fmdiv(m03,mpi,mpi)

      RETURN
    END SUBROUTINE fmpi2
    SUBROUTINE fmprnt(ma)

!  Print MA in base 10 format.

!  FMPRNT can be called directly by the user for easy output
!  in M format.  MA is converted using FMOUT and printed.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC int, log10, max, min, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: j, k, ksave, l, last, lb, nd, nexp
      CHARACTER (20) :: form
! ..
! .. External Subroutines ..
      EXTERNAL fmout
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      namest(ncall) = 'FMPRNT'
      ksave = kflag
      nd = int(real(ndig)*log10(real(mbase))) + 1
      IF (nd<2) nd = 2
      nexp = int(2.0*log10(real(mxbase))) + 6
      lb = max(jform2+nexp,nd+nexp)
      lb = min(lb,lmbuff)
      CALL fmout(ma,cmbuff,lb)
      kflag = ksave
      last = lb + 1
      WRITE (form,90000) kswide - 7
      DO 10 j = 1, lb
        IF (cmbuff(last-j)/=' ' .OR. j==lb) THEN
          l = last - j
          WRITE (kw,form) (cmbuff(k),k=1,l)
          ncall = ncall - 1
          RETURN
        END IF
10    CONTINUE
      ncall = ncall - 1
      RETURN
90000 FORMAT (' (6X,',I3,'A1) ')
    END SUBROUTINE fmprnt
    SUBROUTINE fmpwr(ma,mb,mc)

!  MC = MA ** MB

!  If MB can be expressed exactly as a one word integer, then FMIPWR is
!  used.  This is much faster when MB is small, and using FMIPWR allows
!  MA to be negative.

      IMPLICIT NONE

!             Scratch array usage during FMPWR:   M01 - M06

! .. Intrinsic Functions ..
      INTRINSIC abs, int, log, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: macca, maccb, macmax, mxsave
      INTEGER :: iextra, intmb, j, k, kasave, kfl, kovun, kreslt, kwrnsv, &
        ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmentr, fmeq2, fmexit, fmexp, fmim, fmipwr, fmln, fmmi, &
        fmmpy, fmntr, fmrslt, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
!             Convert MB to an integer before changing NDIG.

      kwrnsv = kwarn
      kwarn = 0
      CALL fmmi(mb,intmb)
      kwarn = kwrnsv
      kfl = kflag

      IF (mblogs/=mbase) CALL fmcons
      IF (abs(ma(1))>mexpab .OR. abs(mb(1))>mexpab .OR. ma(2)<=0) THEN
        CALL fmentr('FMPWR ',ma,mb,2,mc,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        namest(ncall) = 'FMPWR '
        IF (ntrace/=0) CALL fmntr(2,ma,mb,2)
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
        IF (mb(1)==mexpov .OR. mb(1)==mexpun) kovun = 1
        ndsave = ndig
        IF (ncall==1) THEN
          k = max(ngrd52-1,2)
          ndig = max(ndig+k,2)
          IF (ndig>ndg2mx) THEN
            kflag = -9
            CALL fmwarn
            ndig = ndsave
            kreslt = 12
            CALL fmrslt(ma,mb,mc,kreslt)
            IF (ntrace/=0) CALL fmntr(1,mc,mc,1)
            ncall = ncall - 1
            RETURN
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 0
        mxsave = mxexp
        mxexp = mxexp2
      END IF

!             If the exponent is large or the base is very large,
!             raise the precision.

      IF (ma(1)/=0) THEN
        iextra = max(0,int(mb(1))) + int(log(abs(real(ma(1))))/alogmb)
      ELSE
        iextra = max(0,int(mb(1)))
      END IF
      IF (mb(1)-ndig>log(alogmb*real(mxexp2))) THEN
        iextra = 0
      END IF

      ndig = ndig + iextra
      IF (ndig>ndg2mx) THEN
        kflag = -9
        CALL fmwarn
        mc(1) = munkno
        mc(2) = 1
        mc(0) = nint(ndig*alogm2)
        DO 10 j = 2, ndsave
          mc(j+1) = 0
10      CONTINUE
        ndig = ndig - iextra
        CALL fmexit(mc,mc,ndsave,mxsave,kasave,kovun)
        RETURN
      END IF

!             If the exponent is a small integer, call FMIPWR.

      kwrnsv = kwarn
      kwarn = 0

      macca = ma(0)
      maccb = nint(ndig*alogm2)
      CALL fmeq2(ma,m06,ndsave,ndig,0)
      m06(0) = nint(ndig*alogm2)

      IF (kfl==0) THEN
        CALL fmipwr(m06,intmb,mc)
      ELSE IF (m06(2)<=0) THEN
        CALL fmim(0,mc)
        mc(1) = munkno
        mc(2) = 1
        mc(0) = nint(ndig*alogm2)
        kflag = -4
      ELSE
        CALL fmln(m06,m06)
        maccb = mb(0)
        CALL fmeq2(mb,m02,ndsave,ndig,0)
        m02(0) = nint(ndig*alogm2)
        CALL fmmpy(m06,m02,m06)
        CALL fmexp(m06,mc)
      END IF
      kwarn = kwrnsv

!             Round the result and return.

      macmax = nint((ndsave-1)*alogm2+log(real(abs(mc(2))+1))/0.69315)
      mc(0) = min(mc(0),macca,maccb,macmax)
      CALL fmexit(mc,mc,ndsave,mxsave,kasave,kovun)
      RETURN
    END SUBROUTINE fmpwr
    SUBROUTINE fmrdc(ma,mb,jsin,jcos,jswap)

!  Reduce MA using various trigonometric identities to an equivalent
!  angle MB between 0 and 45 degrees.  The reduction is done in radians
!  if KRAD (in common /FMUSER/) is 1, in degrees if KRAD is 0.
!  JSIN and JCOS are returned +1 or -1 and JSWAP is returned to indicate
!  that the sin and cos functions have been interchanged as follows:

!  JSWAP = 0 means   SIN(MA) = JSIN*SIN(MB)
!                    COS(MA) = JCOS*COS(MB)

!  JSWAP = 1 means   SIN(MA) = JSIN*COS(MB)
!                    COS(MA) = JCOS*SIN(MB)

      IMPLICIT NONE

!             Scratch array usage during FMRDC:   M01 - M04

! .. Intrinsic Functions ..
      INTRINSIC abs, int, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: jcos, jsin, jswap
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: x
      INTEGER :: j, kasave, ndsave, ndsv
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmdiv, fmdivi, fmeq, fmeq2, fmi2m, fmint, &
        fmm2dp, fmmpy, fmpi, fmsub, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      jsin = 1
      jcos = 1
      jswap = 0
      ndsave = ndig
      ndig = ndig + max(0,int(ma(1)))

!             If the argument is too big, return UNKNOWN.

      IF (ndig>ndg2mx) THEN
        kflag = -9
        CALL fmwarn
        mb(1) = munkno
        mb(2) = 1
        mb(0) = nint(ndig*alogm2)
        DO 10 j = 2, ndsave
          mb(j+1) = 0
10      CONTINUE
        ndig = ndsave
        RETURN
      END IF
      ma(0) = ma(0) + nint(alogm2*real(max(0,int(ma(1)))))

!             If MA is less than 1/MBASE, no reduction is needed.

      IF (ma(1)<0) THEN
        ndig = ndsave
        CALL fmeq(ma,mb)
        IF (mb(2)<0) THEN
          mb(2) = -mb(2)
          jsin = -1
        END IF
        RETURN
      END IF

      j = 1
      IF (krad==1) THEN
20      IF (mbspi/=mbase .OR. ndigpi<ndig) THEN
          ndsv = ndig
          ndig = min(ndig+2,ndg2mx)
          kasave = kaccsw
          kaccsw = 0
          ncall = ncall + 1
          namest(ncall) = 'NOEQ  '
          CALL fmpi(mpisav)
          ncall = ncall - 1
          kaccsw = kasave
          ndig = ndsv
        END IF
        CALL fmeq2(ma,m04,ndsave,ndig,0)
        IF (ma(2)<0) jsin = -1
        m04(2) = abs(m04(2))
        IF (m04(1)==0) THEN
          CALL fmm2dp(m04,x)
          IF (x<=0.75) THEN
            ndig = ndsave
            CALL fmeq(m04,mb)
            RETURN
          END IF
        END IF
        CALL fmadd(mpisav,mpisav,m02)
        IF (fmcomp(m04,'GE',m02)) THEN
          CALL fmdiv(m04,m02,m01)
          CALL fmint(m01,m01)
          CALL fmmpy(m01,m02,m01)
          CALL fmsub(m04,m01,m04)
        END IF
        CALL fmeq(mpisav,m03)
        IF (fmcomp(m04,'GE',m03)) THEN
          jsin = -jsin
          CALL fmsub(m02,m04,m04)
        END IF
        CALL fmdivi(m02,4,m02)
        IF (fmcomp(m04,'GE',m02)) THEN
          jcos = -jcos
          CALL fmsub(m03,m04,m04)
        END IF
        CALL fmdivi(m03,4,m03)
        IF (fmcomp(m04,'GE',m03)) THEN
          jswap = 1
          CALL fmsub(m02,m04,m04)
        END IF

!             If the reduced argument is close to zero, then
!             cancellation has produced an inaccurate value.
!             Raise NDIG and do the reduction again.

        IF (j==1 .AND. (m04(1)<0 .OR. m04(2)==0)) THEN
          j = 2
          IF (m04(2)==0) THEN
            ndig = min(2*ndig,ndg2mx)
          ELSE
            ndig = ndig - int(m04(1))
          END IF
          IF (ndig>ndg2mx) THEN
            kflag = -9
            CALL fmwarn
            mb(1) = munkno
            mb(2) = 1
            mb(0) = nint(ndig*alogm2)
            DO 30 j = 2, ndsave
              mb(j+1) = 0
30          CONTINUE
            ndig = ndsave
            RETURN
          END IF
          jsin = 1
          jcos = 1
          jswap = 0
          ma(0) = ma(0) + nint(alogm2*real(-m04(1)))
          GO TO 20
        END IF

      ELSE

        CALL fmeq2(ma,m04,ndsave,ndig,0)
        IF (ma(2)<0) jsin = -1
        m04(2) = abs(m04(2))
        IF (m04(1)==0) THEN
          CALL fmm2dp(m04,x)
          IF (x<=44.0) THEN
            ndig = ndsave
            CALL fmeq(m04,mb)
            RETURN
          END IF
        END IF
        CALL fmi2m(360,m02)
        IF (fmcomp(m04,'GE',m02)) THEN
          CALL fmdiv(m04,m02,m01)
          CALL fmint(m01,m01)
          CALL fmmpy(m01,m02,m01)
          CALL fmsub(m04,m01,m04)
        END IF
        CALL fmi2m(180,m03)
        IF (fmcomp(m04,'GE',m03)) THEN
          jsin = -jsin
          CALL fmsub(m02,m04,m04)
        END IF
        CALL fmi2m(90,m02)
        IF (fmcomp(m04,'GE',m02)) THEN
          jcos = -jcos
          CALL fmsub(m03,m04,m04)
        END IF
        CALL fmi2m(45,m03)
        IF (fmcomp(m04,'GE',m03)) THEN
          jswap = 1
          CALL fmsub(m02,m04,m04)
        END IF

      END IF

!             Round the result and return.

      CALL fmeq2(m04,mb,ndig,ndsave,0)
      ndig = ndsave
      RETURN
    END SUBROUTINE fmrdc
    SUBROUTINE fmread(kread,ma)

!  Read MA on unit KREAD.  Multi-line numbers will have '&' as the
!  last nonblank character on all but the last line.  Only one
!  number is allowed on the line(s).

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC max, min, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kread
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: j, lb, ndsave
! ..
! .. Local Arrays ..
      CHARACTER (1) :: line(80)
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmeq2, fminp, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      ncall = ncall + 1
      namest(ncall) = 'FMREAD'
      ndsave = ndig
      ndig = min(ndg2mx,max(ndig+ngrd52,2))
      lb = 0

10    READ (kread,90000,err=30,end=30) line

!             Scan the line and look for '&'

      DO 20 j = 1, 80
        IF (line(j)=='&') GO TO 10
        IF (line(j)/=' ') THEN
          lb = lb + 1
          IF (lb>lmbuff) THEN
            kflag = -8
            GO TO 40
          END IF
          cmbuff(lb) = line(j)
        END IF
20    CONTINUE

      CALL fminp(cmbuff,m01,1,lb)

      CALL fmeq2(m01,ma,ndig,ndsave,0)
      ndig = ndsave
      ncall = ncall - 1
      RETURN

!             If there is an error, return UNKNOWN.

30    kflag = -4
40    CALL fmwarn
      ma(1) = munkno
      ma(2) = 1
      ma(0) = nint(ndig*alogm2)
      DO 50 j = 2, ndig
        ma(j+1) = 0
50    CONTINUE
      ncall = ncall - 1
      RETURN
90000 FORMAT (80A1)
    END SUBROUTINE fmread
    SUBROUTINE fmrnd(mw,nd,nguard,kshift)

!  Round MW to ND digits (base MBASE).

!  MW is non-negative and has ND+NGUARD+KSHIFT digits.

!  NGUARD is the number of guard digits carried.
!  KSHIFT is 1 if a left shift is pending when MW(2)=0.

!  Round to position MW(ND+1+KSHIFT) using the guard digits
!  MW(ND+2+KSHIFT), ..., MW(ND+1+NGUARD+KSHIFT).

!  This routine is designed to be called only from within the FM
!  package.  The user should call FMEQU to round numbers.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC dint, int
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kshift, nd, nguard
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: mw(lmwa)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: m2, mfactr, mkt
      INTEGER :: j, k, kb, l
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (kround==0 .AND. ncall<=1) RETURN
      l = nd + 2 + kshift
      IF (2*(mw(l)+1)<mbase) RETURN
      IF (2*mw(l)>mbase) THEN
        mw(l-1) = mw(l-1) + 1
        mw(l) = 0
        IF (mw(l-1)<mbase) RETURN
        GO TO 40
      END IF

!             If the first guard digit gives a value close to 1/2 then
!             further guard digits must be examined.

      m2 = 2
      IF (int(mbase-dint(mbase/m2)*m2)==0) THEN
        IF (2*mw(l)<mbase) RETURN
        IF (2*mw(l)==mbase) THEN
          IF (nguard>=2) THEN
            IF (mbase>=1000) THEN
              IF (mbase<1000000) THEN
                mfactr = int(0.5D0+0.6883D0*mbase)
              ELSE
                mfactr = int(0.5D0+0.687783D0*mbase)
              END IF
              IF (mw(l+1)==mfactr) RETURN
            END IF
            DO 10 j = 2, nguard
              IF (mw(l+j-1)>0) GO TO 30
10          CONTINUE
          END IF

!                       Round to even.

          IF (int(mw(l-1)-dint(mw(l-1)/m2)*m2)==0) RETURN
        END IF
      ELSE
        IF (2*mw(l)+1==mbase) THEN
          IF (nguard>=2) THEN
            DO 20 j = 2, nguard
              IF (2*(mw(l+j-1)+1)<mbase) RETURN
              IF (2*mw(l+j-1)>mbase) GO TO 30
20          CONTINUE
            RETURN
          END IF
        END IF
      END IF

30    mw(l-1) = mw(l-1) + 1
      mw(l) = 0

!             Check whether there was a carry in the rounded digit.

40    kb = l - 1
      IF (kb>=3) THEN
        k = kb + 1
        DO 50 j = 3, kb
          k = k - 1
          IF (mw(k)<mbase) RETURN
          mkt = dint(mw(k)/mbase)
          mw(k-1) = mw(k-1) + mkt
          mw(k) = mw(k) - mkt*mbase
50      CONTINUE
      END IF

!             If there is a carry in the first digit then the exponent
!             must be adjusted and the number shifted right.

      IF (mw(2)>=mbase) THEN
        IF (kb>=4) THEN
          k = kb + 1
          DO 60 j = 4, kb
            k = k - 1
            mw(k) = mw(k-1)
60        CONTINUE
        END IF

        mkt = dint(mw(2)/mbase)
        IF (kb>=3) mw(3) = mw(2) - mkt*mbase
        mw(2) = mkt
        mw(1) = mw(1) + 1
      END IF

      RETURN
    END SUBROUTINE fmrnd
    SUBROUTINE fmrpwr(ma,ival,jval,mb)

!  MB = MA ** (IVAL/JVAL)   rational exponentiation.

!  This routine is faster than FMPWR when IVAL and JVAL are
!  small integers.

      IMPLICIT NONE

!             Scratch array usage during FMRPWR:   M01 - M03

! .. Intrinsic Functions ..
      INTRINSIC abs, dble, int, log, max, min, mod, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival, jval
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: f, x
      REAL (KIND(0.0D0)) :: ma1, ma2, macca, macmax, mxsave
      REAL :: xval
      INTEGER :: ijsign, invert, ival2, j, jval2, k, kasave, kovun, kreslt, &
        kst, kwrnsv, l, lval, ndsave
! ..
! .. Local Arrays ..
      INTEGER :: nstack(19)
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmdig, fmdiv, fmdivi, fmdpm, fmeq, fmeq2, &
        fmexit, fmgcdi, fmi2m, fmim, fmipwr, fmm2dp, fmmpyi, fmntr, fmntri, &
        fmrslt, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      ncall = ncall + 1
      namest(ncall) = 'FMRPWR'
      IF (ntrace/=0) THEN
        CALL fmntr(2,ma,ma,1)
        CALL fmntri(2,ival,0)
        CALL fmntri(2,jval,0)
      END IF
      kovun = 0
      IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
      ndsave = ndig
      IF (ncall==1) THEN
        xval = max(abs(ival),abs(jval))
        k = int((5.0*real(dlogtn)+2.0*log(xval))/alogmb+2.0)
        ndig = max(ndig+k,2)
      ELSE
        xval = max(abs(ival),abs(jval))
        k = int(log(xval)/alogmb+1.0)
        ndig = ndig + k
      END IF
      IF (ndig>ndg2mx) THEN
        kflag = -9
        CALL fmwarn
        ndig = ndsave
        kreslt = 12
        CALL fmrslt(ma,ma,mb,kreslt)
        IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
        ncall = ncall - 1
        RETURN
      END IF
      kasave = kaccsw
      kaccsw = 0
      mxsave = mxexp
      mxexp = mxexp2

      ma1 = ma(1)
      ma2 = ma(2)
      macca = ma(0)
      CALL fmeq2(ma,m02,ndsave,ndig,0)
      m02(0) = nint(ndig*alogm2)

!             Use GCD-reduced positive exponents.

      ijsign = 1
      ival2 = abs(ival)
      jval2 = abs(jval)
      IF (ival>0 .AND. jval<0) ijsign = -1
      IF (ival<0 .AND. jval>0) ijsign = -1
      IF (ival2>0 .AND. jval2>0) CALL fmgcdi(ival2,jval2)

!             Check for special cases.

10    IF (ma1==munkno .OR. jval2==0 .OR. (ijsign<=0 .AND. ma2==0)) THEN
        CALL fmim(0,mb)
        mb(1) = munkno
        mb(2) = 1
        mb(0) = nint(ndig*alogm2)
        kflag = -4
        GO TO 30
      END IF

      IF (ival2==0) THEN
        CALL fmim(1,mb)
        GO TO 30
      END IF

      IF (jval2==1) THEN
        CALL fmipwr(m02,ijsign*ival2,mb)
        GO TO 30
      END IF

      IF (ma2==0) THEN
        CALL fmeq(ma,mb)
        GO TO 30
      END IF

      IF (ma2<0) THEN
        IF (mod(jval2,2)==0) THEN
          jval2 = 0
          GO TO 10
        END IF
      END IF

      IF (ma1==mexpov) THEN
        IF (ival2<jval2) THEN
          jval2 = 0
          GO TO 10
        END IF
        CALL fmim(0,mb)
        IF (ijsign==1 .AND. ma2>0) THEN
          mb(1) = mexpov
          mb(2) = 1
          mb(0) = nint(ndig*alogm2)
          kflag = -5
        ELSE IF (ijsign==-1 .AND. ma2>0) THEN
          mb(1) = mexpun
          mb(2) = 1
          mb(0) = nint(ndig*alogm2)
          kflag = -6
        ELSE IF (ijsign==1 .AND. ma2<0) THEN
          IF (mod(ival2,2)==0) THEN
            mb(1) = mexpov
            mb(2) = 1
            mb(0) = nint(ndig*alogm2)
            kflag = -5
          ELSE
            mb(1) = mexpov
            mb(2) = -1
            mb(0) = nint(ndig*alogm2)
            kflag = -5
          END IF
        ELSE IF (ijsign==-1 .AND. ma2<0) THEN
          IF (mod(ival2,2)==0) THEN
            mb(1) = mexpun
            mb(2) = 1
            mb(0) = nint(ndig*alogm2)
            kflag = -6
          ELSE
            mb(1) = mexpun
            mb(2) = -1
            mb(0) = nint(ndig*alogm2)
            kflag = -6
          END IF
        END IF
        GO TO 30
      END IF

      IF (ma1==mexpun) THEN
        IF (ival2<jval2) THEN
          jval2 = 0
          GO TO 10
        END IF
        CALL fmim(0,mb)
        IF (ijsign==1 .AND. ma2>0) THEN
          mb(1) = mexpun
          mb(2) = 1
          mb(0) = nint(ndig*alogm2)
          kflag = -6
        ELSE IF (ijsign==-1 .AND. ma2>0) THEN
          mb(1) = mexpov
          mb(2) = 1
          mb(0) = nint(ndig*alogm2)
          kflag = -5
        ELSE IF (ijsign==1 .AND. ma2<0) THEN
          IF (mod(ival2,2)==0) THEN
            mb(1) = mexpun
            mb(2) = 1
            mb(0) = nint(ndig*alogm2)
            kflag = -6
          ELSE
            mb(1) = mexpun
            mb(2) = -1
            mb(0) = nint(ndig*alogm2)
            kflag = -6
          END IF
        ELSE IF (ijsign==-1 .AND. ma2<0) THEN
          IF (mod(ival2,2)==0) THEN
            mb(1) = mexpov
            mb(2) = 1
            mb(0) = nint(ndig*alogm2)
            kflag = -5
          ELSE
            mb(1) = mexpov
            mb(2) = -1
            mb(0) = nint(ndig*alogm2)
            kflag = -5
          END IF
        END IF
        GO TO 30
      END IF

!             Invert MA if MA > 1 and IVAL or JVAL is large.

      invert = 0
      IF (ma(1)>0) THEN
        IF (ival>5 .OR. jval>5) THEN
          invert = 1
          CALL fmi2m(1,m01)
          CALL fmdiv(m01,m02,m02)
        END IF
      END IF

!             Generate the first approximation to ABS(MA)**(1/JVAL2).

      ma1 = m02(1)
      m02(1) = 0
      m02(2) = abs(m02(2))
      CALL fmm2dp(m02,x)
      l = int(ma1/jval2)
      f = ma1/dble(jval2) - l
      x = x**(1.0D0/jval2)*dble(mbase)**f
      CALL fmdpm(x,mb)
      mb(1) = mb(1) + l
      m02(1) = ma1

!             Initialize.

      CALL fmdig(nstack,kst)

!             Newton iteration.

      DO 20 j = 1, kst
        ndig = nstack(j)
        IF (j<kst) ndig = ndig + 1
        lval = jval2 - 1
        CALL fmipwr(mb,lval,m03)
        CALL fmdiv(m02,m03,m03)
        CALL fmmpyi(mb,lval,mb)
        CALL fmadd(mb,m03,mb)
        CALL fmdivi(mb,jval2,mb)
20    CONTINUE

      IF (mb(1)/=munkno .AND. ma2<0) mb(2) = -mb(2)
      CALL fmipwr(mb,ijsign*ival2,mb)
      IF (invert==1) THEN
        CALL fmi2m(1,m01)
        CALL fmdiv(m01,mb,mb)
      END IF

!             Round the result and return.

30    macmax = nint((ndsave-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
      mb(0) = min(macca,macmax)
      kwrnsv = kwarn
      IF (ma1==munkno) kwarn = 0
      CALL fmexit(mb,mb,ndsave,mxsave,kasave,kovun)
      kwarn = kwrnsv
      RETURN
    END SUBROUTINE fmrpwr
    SUBROUTINE fmrslt(ma,mb,mc,kreslt)

!  Handle results that are special cases, such as overflow,
!  underflow, and unknown.

!  MA and MB are the input arguments to an FM subroutine.

!  MC is the result that is returned.

!  KRESLT is the result code from FMARGS.  Result codes handled here:

!   0 - Perform the normal operation
!   1 - The result is the first input argument
!   2 - The result is the second input argument
!   3 - The result is -OVERFLOW
!   4 - The result is +OVERFLOW
!   5 - The result is -UNDERFLOW
!   6 - The result is +UNDERFLOW
!   7 - The result is -1.0
!   8 - The result is +1.0
!  11 - The result is 0.0
!  12 - The result is UNKNOWN

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC min
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kreslt
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccab
      INTEGER :: kfsave
! ..
! .. External Subroutines ..
      EXTERNAL fmeq, fmim
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      kfsave = kflag
      maccab = min(ma(0),mb(0))
      IF (kreslt==1) THEN
        CALL fmeq(ma,mc)
        mc(0) = maccab
        IF (namest(ncall)=='FMADD ' .OR. namest(ncall)=='FMSUB ') THEN
          kflag = 1
        ELSE
          kflag = kfsave
        END IF
        RETURN
      END IF

      IF (kreslt==2) THEN
        CALL fmeq(mb,mc)
        mc(0) = maccab
        IF (namest(ncall)=='FMADD ') THEN
          kflag = 1
        ELSE
          kflag = kfsave
        END IF
        IF (namest(ncall)=='FMSUB ') THEN
          mc(2) = -mc(2)
          kflag = kfsave
        END IF
        RETURN
      END IF

      IF (kreslt==3 .OR. kreslt==4) THEN
        CALL fmim(0,mc)
        mc(1) = mexpov
        mc(2) = 1
        IF (kreslt==3) mc(2) = -1
        mc(0) = maccab
        kflag = kfsave
        RETURN
      END IF

      IF (kreslt==5 .OR. kreslt==6) THEN
        CALL fmim(0,mc)
        mc(1) = mexpun
        mc(2) = 1
        IF (kreslt==5) mc(2) = -1
        mc(0) = maccab
        kflag = kfsave
        RETURN
      END IF

      IF (kreslt==7) THEN
        CALL fmim(-1,mc)
        mc(0) = maccab
        kflag = kfsave
        RETURN
      END IF

      IF (kreslt==8) THEN
        CALL fmim(1,mc)
        mc(0) = maccab
        kflag = kfsave
        RETURN
      END IF

      IF (kreslt==11) THEN
        CALL fmim(0,mc)
        mc(0) = maccab
        kflag = kfsave
        RETURN
      END IF

      IF (kreslt==12 .OR. kreslt<0 .OR. kreslt>15) THEN
        CALL fmim(0,mc)
        mc(1) = munkno
        mc(2) = 1
        mc(0) = maccab
        kflag = kfsave
        RETURN
      END IF

      RETURN
    END SUBROUTINE fmrslt
    SUBROUTINE fmsign(ma,mb,mc)

!  MC = SIGN(MA,MB)

!  MC is set to ABS(MA) if MB is positive or zero,
!     or -ABS(MA) if MB is negative.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: kwrnsv
! ..
! .. External Subroutines ..
      EXTERNAL fmeq, fmim, fmntr
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      kflag = 0
      ncall = ncall + 1
      namest(ncall) = 'FMSIGN'
      IF (ntrace/=0) CALL fmntr(2,ma,mb,2)

      kwrnsv = kwarn
      kwarn = 0
      IF (ma(1)==munkno .OR. mb(1)==munkno) THEN
        CALL fmim(0,mc)
        mc(1) = munkno
        mc(2) = 1
        mc(0) = nint(ndig*alogm2)
        kflag = -4
      ELSE IF (mb(2)>=0) THEN
        CALL fmeq(ma,mc)
        mc(2) = abs(mc(2))
      ELSE
        CALL fmeq(ma,mc)
        mc(2) = -abs(mc(2))
      END IF

      kwarn = kwrnsv
      IF (ntrace/=0) CALL fmntr(1,mc,mc,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmsign
    SUBROUTINE fmsin(ma,mb)

!  MB = SIN(MA)

      IMPLICIT NONE

!             Scratch array usage during FMSIN:   M01 - M04

! .. Intrinsic Functions ..
      INTRINSIC abs, log, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, macca, macmax, mxsave
      INTEGER :: jcos, jsin, jswap, k, kasave, kovun, kreslt, ndsave, ndsv
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmcos2, fmdivi, fmentr, fmeq2, fmexit, fmi2m, fmmpy, &
        fmntr, fmpi, fmrdc, fmrslt, fmsin2, fmsqr, fmsqrt, fmsub, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (abs(ma(1))>mexpab .OR. ma(2)==0) THEN
        CALL fmentr('FMSIN ',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        namest(ncall) = 'FMSIN '
        IF (ntrace/=0) CALL fmntr(2,ma,ma,1)
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
        ndsave = ndig
        IF (ncall==1) THEN
          k = max(ngrd52,2)
          ndig = max(ndig+k,2)
          IF (ndig>ndg2mx) THEN
            kflag = -9
            CALL fmwarn
            ndig = ndsave
            kreslt = 12
            CALL fmrslt(ma,ma,mb,kreslt)
            IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
            ncall = ncall - 1
            RETURN
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 0
        mxsave = mxexp
        mxexp = mxexp2
      END IF

      macca = ma(0)
      ma2 = ma(2)
      CALL fmeq2(ma,mb,ndsave,ndig,0)
      mb(0) = nint(ndig*alogm2)
      mb(2) = abs(mb(2))

!             Reduce the argument, convert to radians if the input is
!             in degrees, and evaluate the function.

      CALL fmrdc(mb,mb,jsin,jcos,jswap)
      IF (mb(1)==munkno) GO TO 10
      IF (krad==0) THEN
        IF (mbspi/=mbase .OR. ndigpi<ndig) THEN
          ndsv = ndig
          ndig = min(ndig+2,ndg2mx)
          ncall = ncall + 1
          namest(ncall) = 'NOEQ  '
          CALL fmpi(mpisav)
          ncall = ncall - 1
          ndig = ndsv
        END IF
        CALL fmmpy(mb,mpisav,mb)
        CALL fmdivi(mb,180,mb)
      END IF
      IF (mb(1)/=munkno) THEN
        IF (jswap==0) THEN
          IF (mb(1)<0 .OR. ndig<=50) THEN
            CALL fmsin2(mb,mb)
          ELSE
            CALL fmcos2(mb,mb)
            CALL fmi2m(1,m03)
            CALL fmsqr(mb,mb)
            CALL fmsub(m03,mb,mb)
            CALL fmsqrt(mb,mb)
          END IF
        ELSE
          CALL fmcos2(mb,mb)
        END IF
      END IF

!             Append the sign, round, and return.

      IF (jsin==-1 .AND. mb(1)/=munkno) mb(2) = -mb(2)
10    macmax = nint((ndsave-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
      mb(0) = min(mb(0),macca,macmax)
      IF (ma2<0 .AND. mb(1)/=munkno) mb(2) = -mb(2)
      CALL fmexit(mb,mb,ndsave,mxsave,kasave,kovun)
      RETURN
    END SUBROUTINE fmsin
    SUBROUTINE fmsin2(ma,mb)

!  Internal subroutine for MB = SIN(MA) where 0.LE.MA.LE.1.

      IMPLICIT NONE

!             Scratch array usage during FMSIN2:   M01 - M04

!             LJSUMS = 8*(LUNPCK+1) allows for up to eight concurrent
!             sums.  Increasing this value will begin to improve the
!             speed of SIN when the base is large and precision exceeds
!             about 1,500 decimal digits.

! .. Intrinsic Functions ..
      INTRINSIC int, log, max, min, nint, real, sqrt
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL :: alog3, alogt, b, t, tj
      REAL (KIND(0.0D0)) :: maxval
      INTEGER :: j, j2, k, k2, kpt, kthree, kwrnsv, l, l2, large, n2, nbot, &
        ndsav1, ndsave, nterm
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmdivi, fmeq, fmeq2, fmi2m, fmipwr, fmmpy, &
        fmmpyi, fmsqr, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mjsums(0:ljsums), &
        mlbsav(0:lunpck), mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), &
        mln4(0:lunpck), mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmsums/mjsums
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (ma(2)==0) THEN
        CALL fmeq(ma,mb)
        RETURN
      END IF
      ndsave = ndig
      kwrnsv = kwarn
      kwarn = 0

!             Use the direct series
!                  SIN(X) = X - X**3/3! + X**5/5! - ...

!             The argument will be divided by 3**K2 before the series
!             is summed.  The series will be added as J2 concurrent
!             series.  The approximately optimal values of K2 and J2
!             are now computed to try to minimize the time required.
!             N2/2 is the approximate number of terms of the series
!             that will be needed, and L2 guard digits will be carried.

      b = real(mbase)
      k = ngrd52
      t = max(ndig-k,2)
      alog3 = log(3.0)
      alogt = log(t)
      tj = 0.05*alogmb*t**0.3333 + 1.85
      j2 = int(tj)
      j2 = max(1,min(j2,ljsums/ndg2mx))
      k2 = int(0.1*sqrt(t*alogmb/tj)-0.05*alogt+2.5)

      l = int(-(real(ma(1))*alogmb+log(real(ma(2))/b+ &
        real(ma(3))/(b*b)))/alog3-0.3)
      k2 = k2 - l
      IF (l<0) l = 0
      IF (k2<0) THEN
        k2 = 0
        j2 = int(.43*sqrt(t*alogmb/(alogt+real(l)*alog3))+.33)
      END IF
      IF (j2<=1) j2 = 1

      n2 = int(t*alogmb/(alogt+real(l)*alog3))
      l2 = int(log(real(n2)+3.0**k2)/alogmb)
      ndig = ndig + l2
      IF (ndig>ndg2mx) THEN
        kflag = -9
        CALL fmwarn
        mb(1) = munkno
        mb(2) = 1
        mb(0) = nint(ndig*alogm2)
        DO 10 j = 2, ndsave
          mb(j+1) = 0
10      CONTINUE
        ndig = ndsave
        kwarn = kwrnsv
        RETURN
      END IF
      ndsav1 = ndig

!             Divide the argument by 3**K2.

      CALL fmeq2(ma,m02,ndsave,ndig,0)
      kthree = 1
      maxval = mxbase/3
      IF (k2>0) THEN
        DO 20 j = 1, k2
          kthree = 3*kthree
          IF (kthree>maxval) THEN
            CALL fmdivi(m02,kthree,m02)
            kthree = 1
          END IF
20      CONTINUE
        IF (kthree>1) CALL fmdivi(m02,kthree,m02)
      END IF

!             Split into J2 concurrent sums and reduce NDIG while
!             computing each term in the sum as the terms get smaller.

      CALL fmeq(m02,m03)
      nterm = 1
      DO 30 j = 1, j2
        nbot = nterm*(nterm-1)
        IF (nbot>1) CALL fmdivi(m03,nbot,m03)
        nterm = nterm + 2
        kpt = (j-1)*(ndig+2)
        CALL fmeq(m03,mjsums(kpt))
        m03(2) = -m03(2)
30    CONTINUE
      CALL fmsqr(m02,m02)
      IF (m02(1)<-ndig) GO TO 60
      CALL fmipwr(m02,j2,mb)

40    CALL fmmpy(m03,mb,m03)
      large = int(intmax/nterm)
      DO 50 j = 1, j2
        nbot = nterm*(nterm-1)
        IF (nterm>large .OR. nbot>mxbase) THEN
          CALL fmdivi(m03,nterm,m03)
          nbot = nterm - 1
          CALL fmdivi(m03,nbot,m03)
        ELSE
          CALL fmdivi(m03,nbot,m03)
        END IF
        kpt = (j-1)*(ndsav1+2)
        ndig = ndsav1
        CALL fmadd(mjsums(kpt),m03,mjsums(kpt))
        IF (kflag/=0) GO TO 60
        ndig = ndsav1 - int(mjsums(kpt+1)-m03(1))
        IF (ndig<2) ndig = 2
        m03(2) = -m03(2)
        nterm = nterm + 2
50    CONTINUE
      GO TO 40

!             Next put the J2 separate sums back together.

60    kflag = 0
      kpt = (j2-1)*(ndig+2)
      CALL fmeq(mjsums(kpt),mb)
      IF (j2>=2) THEN
        DO 70 j = 2, j2
          CALL fmmpy(m02,mb,mb)
          kpt = (j2-j)*(ndig+2)
          CALL fmadd(mb,mjsums(kpt),mb)
70      CONTINUE
      END IF

!             Reverse the effect of reducing the argument to
!             compute SIN(MA).

      ndig = ndsav1
      IF (k2>0) THEN
        CALL fmi2m(3,m02)
        DO 80 j = 1, k2
          CALL fmsqr(mb,m03)
          CALL fmmpyi(m03,-4,m03)
          CALL fmadd(m02,m03,m03)
          CALL fmmpy(m03,mb,mb)
80      CONTINUE
      END IF

      CALL fmeq2(mb,mb,ndsav1,ndsave,1)
      ndig = ndsave
      kwarn = kwrnsv

      RETURN
    END SUBROUTINE fmsin2
    SUBROUTINE fmsinh(ma,mb)

!  MB = SINH(MA)

      IMPLICIT NONE

!             Scratch array usage during FMSINH:   M01 - M03

! .. Intrinsic Functions ..
      INTRINSIC abs, int, log, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, macca, macmax, mxsave
      INTEGER :: k, kasave, kovun, kreslt, ndsave, nmethd
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmcsh2, fmdiv, fmdivi, fmentr, fmeq, fmeq2, fmexit, &
        fmexp, fmi2m, fmntr, fmrslt, fmsnh2, fmsqr, fmsqrt, fmsub, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (abs(ma(1))>mexpab) THEN
        CALL fmentr('FMSINH',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        namest(ncall) = 'FMSINH'
        IF (ntrace/=0) CALL fmntr(2,ma,ma,1)
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
        ndsave = ndig
        IF (ncall==1) THEN
          k = max(ngrd52,2)
          ndig = max(ndig+k,2)
          IF (ndig>ndg2mx) THEN
            kflag = -9
            CALL fmwarn
            ndig = ndsave
            kreslt = 12
            CALL fmrslt(ma,ma,mb,kreslt)
            IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
            ncall = ncall - 1
            RETURN
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 0
        mxsave = mxexp
        mxexp = mxexp2
      END IF

      macca = ma(0)
      ma2 = ma(2)
      CALL fmeq2(ma,mb,ndsave,ndig,0)
      mb(0) = nint(ndig*alogm2)
      mb(2) = abs(mb(2))
      IF (ma2==0) THEN
        CALL fmeq(ma,mb)
        GO TO 20
      END IF

!             Use a series for small arguments, FMEXP for large ones.

      IF (mb(1)==munkno) GO TO 20
      IF (mbase>99) THEN
        IF (mb(1)<=0) THEN
          nmethd = 1
        ELSE IF (mb(1)>=2) THEN
          nmethd = 2
        ELSE IF (abs(mb(2))<10) THEN
          nmethd = 1
        ELSE
          nmethd = 2
        END IF
      ELSE
        IF (mb(1)<=0) THEN
          nmethd = 1
        ELSE
          nmethd = 2
        END IF
      END IF

      IF (nmethd==2) GO TO 10
      IF (mb(1)<0 .OR. ndig<=50) THEN
        CALL fmsnh2(mb,mb)
      ELSE
        CALL fmcsh2(mb,mb)
        CALL fmi2m(1,m03)
        CALL fmsqr(mb,mb)
        CALL fmsub(mb,m03,mb)
        CALL fmsqrt(mb,mb)
      END IF
      GO TO 20

10    CALL fmexp(mb,mb)
      IF (mb(1)==mexpov) THEN
        GO TO 20
      ELSE IF (mb(1)==mexpun) THEN
        mb(1) = mexpov
        GO TO 20
      END IF
      IF (int(mb(1))<=(ndig+1)/2) THEN
        CALL fmi2m(1,m01)
        CALL fmdiv(m01,mb,m01)
        CALL fmsub(mb,m01,mb)
      END IF
      CALL fmdivi(mb,2,mb)

!             Round and return.

20    macmax = nint((ndsave-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
      mb(0) = min(mb(0),macca,macmax)
      IF (ma2<0 .AND. mb(1)/=munkno) mb(2) = -mb(2)
      CALL fmexit(mb,mb,ndsave,mxsave,kasave,kovun)
      RETURN
    END SUBROUTINE fmsinh
    SUBROUTINE fmsnh2(ma,mb)

!  Internal subroutine for MB = SINH(MA).

      IMPLICIT NONE

!             Scratch array usage during FMSNH2:   M01 - M03

!             LJSUMS = 8*(LUNPCK+1) allows for up to eight concurrent
!             sums.  Increasing this value will begin to improve the
!             speed of SINH when the base is large and precision exceeds
!             about 1,500 decimal digits.

! .. Intrinsic Functions ..
      INTRINSIC int, log, max, min, nint, real, sqrt
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL :: alog3, alogt, b, t, tj
      REAL (KIND(0.0D0)) :: maxval
      INTEGER :: j, j2, k, k2, kpt, kthree, kwrnsv, l, l2, large, n2, nbot, &
        ndsav1, ndsave, nterm
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmdivi, fmeq, fmeq2, fmi2m, fmipwr, fmmpy, &
        fmmpyi, fmsqr, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mjsums(0:ljsums), &
        mlbsav(0:lunpck), mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), &
        mln4(0:lunpck), mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmsums/mjsums
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (ma(2)==0) THEN
        CALL fmeq(ma,mb)
        RETURN
      END IF
      ndsave = ndig
      kwrnsv = kwarn
      kwarn = 0

!             Use the direct series
!                  SINH(X) = X + X**3/3! + X**5/5! - ...

!             The argument will be divided by 3**K2 before the series
!             is summed.  The series will be added as J2 concurrent
!             series.  The approximately optimal values of K2 and J2
!             are now computed to try to minimize the time required.
!             N2/2 is the approximate number of terms of the series
!             that will be needed, and L2 guard digits will be carried.

      b = real(mbase)
      k = ngrd52
      t = max(ndig-k,2)
      alog3 = log(3.0)
      alogt = log(t)
      tj = 0.05*alogmb*t**0.3333 + 1.85
      j2 = int(tj)
      j2 = max(1,min(j2,ljsums/ndg2mx))
      k2 = int(0.1*sqrt(t*alogmb/tj)-0.05*alogt+2.5)

      l = int(-(real(ma(1))*alogmb+log(real(ma(2))/b+ &
        real(ma(3))/(b*b)))/alog3-0.3)
      k2 = k2 - l
      IF (l<0) l = 0
      IF (k2<0) THEN
        k2 = 0
        j2 = int(.43*sqrt(t*alogmb/(alogt+real(l)*alog3))+.33)
      END IF
      IF (j2<=1) j2 = 1

      n2 = int(t*alogmb/(alogt+real(l)*alog3))
      l2 = int(log(real(n2)+3.0**k2)/alogmb)
      ndig = ndig + l2
      IF (ndig>ndg2mx) THEN
        kflag = -9
        CALL fmwarn
        mb(1) = munkno
        mb(2) = 1
        mb(0) = nint(ndig*alogm2)
        DO 10 j = 2, ndsave
          mb(j+1) = 0
10      CONTINUE
        ndig = ndsave
        kwarn = kwrnsv
        RETURN
      END IF
      ndsav1 = ndig

!             Divide the argument by 3**K2.

      CALL fmeq2(ma,m02,ndsave,ndig,0)
      kthree = 1
      maxval = mxbase/3
      IF (k2>0) THEN
        DO 20 j = 1, k2
          kthree = 3*kthree
          IF (kthree>maxval) THEN
            CALL fmdivi(m02,kthree,m02)
            kthree = 1
          END IF
20      CONTINUE
        IF (kthree>1) CALL fmdivi(m02,kthree,m02)
      END IF

!             Split into J2 concurrent sums and reduce NDIG while
!             computing each term in the sum as the terms get smaller.

      CALL fmeq(m02,m03)
      nterm = 1
      DO 30 j = 1, j2
        nbot = nterm*(nterm-1)
        IF (nbot>1) CALL fmdivi(m03,nbot,m03)
        nterm = nterm + 2
        kpt = (j-1)*(ndig+2)
        CALL fmeq(m03,mjsums(kpt))
30    CONTINUE
      CALL fmsqr(m02,m02)
      IF (m02(1)<-ndig) GO TO 60
      CALL fmipwr(m02,j2,mb)

40    CALL fmmpy(m03,mb,m03)
      large = int(intmax/nterm)
      DO 50 j = 1, j2
        nbot = nterm*(nterm-1)
        IF (nterm>large .OR. nbot>mxbase) THEN
          CALL fmdivi(m03,nterm,m03)
          nbot = nterm - 1
          CALL fmdivi(m03,nbot,m03)
        ELSE
          CALL fmdivi(m03,nbot,m03)
        END IF
        kpt = (j-1)*(ndsav1+2)
        ndig = ndsav1
        CALL fmadd(mjsums(kpt),m03,mjsums(kpt))
        IF (kflag/=0) GO TO 60
        ndig = ndsav1 - int(mjsums(kpt+1)-m03(1))
        IF (ndig<2) ndig = 2
        nterm = nterm + 2
50    CONTINUE
      GO TO 40

!             Next put the J2 separate sums back together.

60    kflag = 0
      kpt = (j2-1)*(ndig+2)
      CALL fmeq(mjsums(kpt),mb)
      IF (j2>=2) THEN
        DO 70 j = 2, j2
          CALL fmmpy(m02,mb,mb)
          kpt = (j2-j)*(ndig+2)
          CALL fmadd(mb,mjsums(kpt),mb)
70      CONTINUE
      END IF

!             Reverse the effect of reducing the argument to
!             compute SINH(MA).

      ndig = ndsav1
      IF (k2>0) THEN
        CALL fmi2m(3,m02)
        DO 80 j = 1, k2
          CALL fmsqr(mb,m03)
          CALL fmmpyi(m03,4,m03)
          CALL fmadd(m02,m03,m03)
          CALL fmmpy(m03,mb,mb)
80      CONTINUE
      END IF

      CALL fmeq2(mb,mb,ndsav1,ndsave,1)
      ndig = ndsave
      kwarn = kwrnsv

      RETURN
    END SUBROUTINE fmsnh2
    SUBROUTINE fmsp2m(x,ma)

!  MA = X

!  Convert a single precision number to FM format.

!  In general the relative accuracy of the number returned is only
!  the relative accuracy of a machine precision number.  This may be
!  true even if X can be represented exactly in the machine floating
!  point number system.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dble, int
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL :: x
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: xdp, y, yt
      INTEGER :: k
! ..
! .. External Subroutines ..
      EXTERNAL fmdivi, fmdm, fmim, fmntr, fmntrr
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      namest(ncall) = 'FMSP2M'
      xdp = dble(x)
      IF (ntrace/=0) CALL fmntrr(2,xdp,1)

!             Check to see if X is exactly a small integer.  If so,
!             converting as an integer is better.
!             Also see if X is exactly a small integer divided by
!             a small power of two.

      y = mxexp2
      IF (abs(xdp)<y) THEN
        k = int(xdp)
        y = k
        IF (y==xdp) THEN
          CALL fmim(k,ma)
          GO TO 10
        END IF
      END IF
      IF (abs(xdp)<1.0D0) THEN
        y = 4096.0D0*xdp
        k = int(y)
        yt = k
        IF (y==yt) THEN
          CALL fmim(k,ma)
          CALL fmdivi(ma,4096,ma)
          GO TO 10
        END IF
      END IF

      CALL fmdm(xdp,ma)

10    IF (ntrace/=0) CALL fmntr(1,ma,ma,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmsp2m
    SUBROUTINE fmsqr(ma,mb)

!  MB = MA*MA    Faster than using FMMPY.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. External Subroutines ..
      EXTERNAL fmntr, fmsqr2
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (ntrace/=0) THEN
        namest(ncall) = 'FMSQR '
        CALL fmntr(2,ma,ma,1)

        CALL fmsqr2(ma,mb)

        CALL fmntr(1,mb,mb,1)
      ELSE
        CALL fmsqr2(ma,mb)
      END IF
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmsqr
    SUBROUTINE fmsqr2(ma,mb)

!  MB = MA*MA.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dint, int, log, max, min, mod, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, macca, maxmax, maxmwa, mbj, mbkj, mbm1, mbnorm, md2b, &
        mk, mka, mkt, mmax, mr, mt
      INTEGER :: j, jm1, k, kb, ki, kj, kl, knz, kovun, kshift, kwa, l, n1, &
        nguard
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmeq, fmmove, fmmpy2, fmrnd, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (abs(ma(1))>mexpab .OR. kdebug==1) THEN
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
        IF (ma(1)==munkno) kovun = 2
        ncall = ncall + 1
        CALL fmmpy2(ma,ma,mb)
        ncall = ncall - 1
        IF ((kflag<0 .AND. kovun==0) .OR. (kflag==-4 .AND. kovun==1)) THEN
          namest(ncall) = 'FMSQR '
          CALL fmwarn
        END IF
        GO TO 140
      ELSE IF (ma(2)==0) THEN
        CALL fmeq(ma,mb)
        GO TO 140
      END IF
      kflag = 0
      maxmax = 0

      macca = ma(0)
      ma2 = ma(2)
      n1 = ndig + 1
      mwa(1) = ma(1) + ma(1)

!             NGUARD is the number of guard digits used.

      IF (ncall>1) THEN
        nguard = ngrd22
        IF (nguard>ndig) nguard = ndig
      ELSE
        nguard = ngrd52
        IF (nguard>ndig) nguard = ndig
      END IF
      IF (ma(2)*ma(2)<mbase .AND. nguard<3) nguard = 3

      l = n1 + nguard
      mwa(l+1) = 0
      ma(2) = abs(ma(2))

!             The multiplication loop begins here.

!             MBNORM is the minimum number of digits that can be
!                    multiplied before normalization is required.
!             MAXMWA is an upper bound on the size of values in MWA
!                    divided by (MBASE-1).  It is used to determine
!                    whether to normalize before the next digit is
!                    multiplied.

      mbm1 = mbase - 1
      mbnorm = dint(maxint/(mbm1*mbm1))
      mmax = intmax - mbase
      mmax = min(dint(maxint/mbm1-mbm1),mmax)
      IF (mbnorm>1) THEN
        mbj = ma(2)

!             Count the trailing zeros in MA.

        IF (ma(n1)/=0) THEN
          knz = n1
        ELSE
          DO 10 j = ndig, 2, -1
            IF (ma(j)/=0) THEN
              knz = j
              GO TO 20
            END IF
10        CONTINUE
        END IF

20      mwa(2) = 0
        mwa(3) = 0
        DO 30 k = ndig + 2, l
          mwa(k) = 0
30      CONTINUE

!             (Inner Loop)

        DO 40 k = 3, n1
          mwa(k+1) = ma(k)*mbj
40      CONTINUE
        maxmwa = mbj
        DO 70 j = 3, l/2
          mbj = ma(j)
          IF (mbj/=0) THEN
            maxmwa = maxmwa + mbj
            jm1 = j - 1
            kl = min(knz,l-jm1)

!                       Major (Inner Loop)

            DO 50 k = 2*j, jm1 + kl
              mwa(k) = mwa(k) + ma(k-jm1)*mbj
50          CONTINUE
          END IF

          IF (maxmwa>mmax) THEN
            maxmax = max(maxmax,maxmwa)
            maxmwa = 0

!                       Normalization is only required for the
!                       range of digits currently changing in MWA.

            DO 60 kb = jm1 + kl, 2*j, -1
              mkt = int(mwa(kb)/mbase)
              mwa(kb-1) = mwa(kb-1) + mkt
              mwa(kb) = mwa(kb) - mkt*mbase
60          CONTINUE
          END IF
70      CONTINUE

!             Double MWA, add the square terms, and perform
!             the final normalization.  (Inner Loop)

        IF (2*max(maxmax,maxmwa)+mbase>mmax) THEN
          DO 80 kb = l, 4, -1
            mkt = int(mwa(kb)/mbase)
            mwa(kb-1) = mwa(kb-1) + mkt
            mwa(kb) = mwa(kb) - mkt*mbase
80        CONTINUE
        END IF

        DO 90 j = 3, l - 1, 2
          mka = ma((j+1)/2)
          mwa(j) = 2*mwa(j) + mka*mka
          mwa(j+1) = 2*mwa(j+1)
90      CONTINUE
        IF (mod(l,2)==1) THEN
          mka = ma((l+1)/2)
          mwa(l) = 2*mwa(l) + mka*mka
        END IF

        DO 100 kb = l, 3, -1
          mkt = int(mwa(kb)/mbase)
          mwa(kb-1) = mwa(kb-1) + mkt
          mwa(kb) = mwa(kb) - mkt*mbase
100     CONTINUE

      ELSE

!             If normalization must be done for each digit, combine
!             the two loops and normalize as the digits are multiplied.

        DO 110 j = 2, l
          mwa(j) = 0
110     CONTINUE
        kj = ndig + 2
        DO 130 j = 2, n1
          kj = kj - 1
          mbkj = ma(kj)
          IF (mbkj==0) GO TO 130
          kl = l - kj + 1
          IF (kl>n1) kl = n1
          ki = kl + 2
          kwa = kl + kj + 1
          mk = 0
          DO 120 k = 2, kl
            mt = ma(ki-k)*mbkj + mwa(kwa-k) + mk
            mk = int(mt/mbase)
            mwa(kwa-k) = mt - mbase*mk
120       CONTINUE
          mwa(kwa-kl-1) = mk
130     CONTINUE

      END IF

!             Set KSHIFT = 1 if a shift left is necessary.

      IF (mwa(2)==0) THEN
        kshift = 1
      ELSE
        kshift = 0
      END IF

!             The multiplication is complete.
!             Round the result and move it to MB.

      ma(2) = ma2
      mr = 2*mwa(ndig+2+kshift) + 1
      IF (mr>=mbase) THEN
        IF (mr-1>mbase .AND. mwa(n1+kshift)<mbase-1) THEN
          IF (kround/=0 .OR. ncall>1) THEN
            mwa(n1+kshift) = mwa(n1+kshift) + 1
            mwa(n1+1+kshift) = 0
          END IF
        ELSE
          CALL fmrnd(mwa,ndig,nguard,kshift)
        END IF
      END IF
      CALL fmmove(mwa,mb)

      IF (kflag<0) THEN
        namest(ncall) = 'FMSQR '
        CALL fmwarn
      END IF

      IF (kaccsw==1) THEN
        md2b = nint((ndig-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
        mb(0) = min(macca,md2b)
      ELSE
        mb(0) = macca
      END IF
140   RETURN
    END SUBROUTINE fmsqr2
    SUBROUTINE fmsqrt(ma,mb)

!  MB = SQRT(MA)

      IMPLICIT NONE

!             Scratch array usage during FMSQRT:   M01 - M02

! .. Intrinsic Functions ..
      INTRINSIC abs, int, log, max, min, mod, nint, real, sqrt
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma1, macca, md2b, mke, mxsave
      REAL (KIND(0.0D0)) :: x, xb
      INTEGER :: j, k, kasave, kma1, kovun, kreslt, kst, ndsave
! ..
! .. Local Arrays ..
      INTEGER :: nstack(19)
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmdig, fmdiv, fmdivi, fmdpm, fmentr, fmeq2, &
        fmexit, fmm2dp, fmntr, fmrslt, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (abs(ma(1))>mexpab .OR. ma(2)<=0) THEN
        CALL fmentr('FMSQRT',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        IF (ntrace/=0) THEN
          namest(ncall) = 'FMSQRT'
          CALL fmntr(2,ma,ma,1)
        END IF
        ndsave = ndig
        IF (ncall==1) THEN
          k = max(ngrd52-1,2)
          ndig = max(ndig+k,2)
          IF (ndig>ndg2mx) THEN
            namest(ncall) = 'FMSQRT'
            kflag = -9
            CALL fmwarn
            ndig = ndsave
            kreslt = 12
            CALL fmrslt(ma,ma,mb,kreslt)
            IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
            ncall = ncall - 1
            RETURN
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 0
        mxsave = mxexp
        mxexp = mxexp2
      END IF

      ma1 = ma(1)

      macca = ma(0)
      CALL fmeq2(ma,m02,ndsave,ndig,0)
      m02(0) = nint(ndig*alogm2)

!             Generate the first approximation.

      m02(1) = 0
      CALL fmm2dp(m02,x)
      x = sqrt(x)
      mke = ma1/2
      kma1 = int(abs(ma1))
      IF (mod(kma1,2)==1) THEN
        xb = mbase
        x = x*sqrt(xb)
        mke = (ma1-1)/2
      END IF
      CALL fmdpm(x,mb)
      mb(1) = mb(1) + mke

!             Initialize.

      m02(1) = ma1
      CALL fmdig(nstack,kst)

!             Newton iteration.

      DO 10 j = 1, kst
        ndig = nstack(j)
        CALL fmdiv(m02,mb,m01)
        CALL fmadd(mb,m01,mb)
        CALL fmdivi(mb,2,mb)
10    CONTINUE

!             Round the result and return.

      IF (kasave==1) THEN
        md2b = nint((ndsave-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
        mb(0) = min(macca,md2b)
      ELSE
        mb(0) = macca
      END IF
      CALL fmexit(mb,mb,ndsave,mxsave,kasave,0)
      RETURN
    END SUBROUTINE fmsqrt
    SUBROUTINE fmst2m(string,ma)

!  MA = STRING

!  Convert a character string to FM format.
!  This is often more convenient than using FMINP, which converts an
!  array of CHARACTER*1 values.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC len
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (*) :: string
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: j, lb
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fminp
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      ncall = ncall + 1
      namest(ncall) = 'FMST2M'
      lb = len(string)

      DO 10 j = 1, lb
        cmbuff(j) = string(j:j)
10    CONTINUE

      CALL fminp(cmbuff,ma,1,lb)

      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmst2m
    SUBROUTINE fmsub(ma,mb,mc)

!  MC = MA - MB

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: kflg1
! ..
! .. External Subroutines ..
      EXTERNAL fmadd2, fmntr
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (ntrace/=0) THEN
        namest(ncall) = 'FMSUB '
        IF (ntrace/=0) CALL fmntr(2,ma,mb,2)

        kflg1 = 0
        IF (mb(1)>ma(1) .OR. ma(2)==0) kflg1 = 1
        IF (mb(2)==0) kflg1 = 0

!             FMADD2 will negate MB and add.

        ksub = 1
        CALL fmadd2(ma,mb,mc)
        ksub = 0

!             If MA was smaller than MB, then KFLAG = 1 returned from
!             FMADD means the result from FMSUB is the opposite of the
!             input argument of larger magnitude, so reset KFLAG.

        IF (kflag==1 .AND. kflg1==1) kflag = 0

        IF (ntrace/=0) CALL fmntr(1,mc,mc,1)
      ELSE
        kflg1 = 0
        IF (mb(1)>ma(1) .OR. ma(2)==0) kflg1 = 1
        IF (mb(2)==0) kflg1 = 0
        ksub = 1
        CALL fmadd2(ma,mb,mc)
        ksub = 0
        IF (kflag==1 .AND. kflg1==1) kflag = 0
      END IF
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmsub
    SUBROUTINE fmtan(ma,mb)

!  MB = TAN(MA)

      IMPLICIT NONE

!             Scratch array usage during FMTAN:   M01 - M04

! .. Intrinsic Functions ..
      INTRINSIC abs, log, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, macca, macmax, mxsave
      INTEGER :: jcos, jsin, jswap, k, kasave, kovun, kreslt, ndsave, ndsv
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmcos2, fmdiv, fmdivi, fmentr, fmeq2, fmexit, fmi2m, &
        fmim, fmmpy, fmntr, fmpi, fmrdc, fmrslt, fmsin2, fmsqr, fmsqrt, fmsub, &
        fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (abs(ma(1))>mexpab .OR. ma(2)==0) THEN
        CALL fmentr('FMTAN ',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        namest(ncall) = 'FMTAN '
        IF (ntrace/=0) CALL fmntr(2,ma,ma,1)
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
        ndsave = ndig
        IF (ncall==1) THEN
          k = max(ngrd52-1,2)
          ndig = max(ndig+k,2)
          IF (ndig>ndg2mx) THEN
            kflag = -9
            CALL fmwarn
            ndig = ndsave
            kreslt = 12
            CALL fmrslt(ma,ma,mb,kreslt)
            IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
            ncall = ncall - 1
            RETURN
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 0
        mxsave = mxexp
        mxexp = mxexp2
      END IF

      macca = ma(0)
      ma2 = ma(2)
      CALL fmeq2(ma,mb,ndsave,ndig,0)
      mb(0) = nint(ndig*alogm2)
      mb(2) = abs(mb(2))

!             Reduce the argument, convert to radians if the input is
!             in degrees, and evaluate the function.

      CALL fmrdc(mb,mb,jsin,jcos,jswap)
      IF (mb(1)==munkno) GO TO 10
      IF (mb(2)==0) THEN
        IF (jswap==1) THEN
          CALL fmim(0,mb)
          mb(1) = munkno
          mb(2) = 1
          mb(0) = nint(ndig*alogm2)
          kflag = -4
          CALL fmwarn
        END IF
        GO TO 10
      END IF
      IF (krad==0) THEN
        IF (mbspi/=mbase .OR. ndigpi<ndig) THEN
          ndsv = ndig
          ndig = min(ndig+2,ndg2mx)
          ncall = ncall + 1
          namest(ncall) = 'NOEQ  '
          CALL fmpi(mpisav)
          ncall = ncall - 1
          ndig = ndsv
        END IF
        CALL fmmpy(mb,mpisav,mb)
        CALL fmdivi(mb,180,mb)
      END IF
      IF (mb(1)/=munkno) THEN
        IF (jswap==0) THEN
          IF (mb(1)<0) THEN
            CALL fmsin2(mb,mb)
            mb(2) = jsin*mb(2)
            CALL fmsqr(mb,m03)
            CALL fmi2m(1,m02)
            CALL fmsub(m02,m03,m03)
            CALL fmsqrt(m03,m04)
            m04(2) = jcos*m04(2)
            CALL fmdiv(mb,m04,mb)
          ELSE
            CALL fmcos2(mb,mb)
            mb(2) = jcos*mb(2)
            CALL fmsqr(mb,m03)
            CALL fmi2m(1,m02)
            CALL fmsub(m02,m03,m03)
            CALL fmsqrt(m03,m04)
            m04(2) = jsin*m04(2)
            CALL fmdiv(m04,mb,mb)
          END IF
        ELSE
          IF (mb(1)<0) THEN
            CALL fmsin2(mb,mb)
            mb(2) = jcos*mb(2)
            CALL fmsqr(mb,m03)
            CALL fmi2m(1,m02)
            CALL fmsub(m02,m03,m03)
            CALL fmsqrt(m03,m04)
            m04(2) = jsin*m04(2)
            CALL fmdiv(m04,mb,mb)
          ELSE
            CALL fmcos2(mb,mb)
            mb(2) = jsin*mb(2)
            CALL fmsqr(mb,m03)
            CALL fmi2m(1,m02)
            CALL fmsub(m02,m03,m03)
            CALL fmsqrt(m03,m04)
            m04(2) = jcos*m04(2)
            CALL fmdiv(mb,m04,mb)
          END IF
        END IF
      END IF

!             Round and return.

10    macmax = nint((ndsave-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
      mb(0) = min(mb(0),macca,macmax)
      IF (ma2<0 .AND. mb(1)/=munkno) mb(2) = -mb(2)
      CALL fmexit(mb,mb,ndsave,mxsave,kasave,kovun)
      RETURN
    END SUBROUTINE fmtan
    SUBROUTINE fmtanh(ma,mb)

!  MB = TANH(MA)

      IMPLICIT NONE

!             Scratch array usage during FMTANH:   M01 - M03

! .. Intrinsic Functions ..
      INTRINSIC abs, int, log, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, macca, macmax, mxsave
      REAL :: x, xt
      INTEGER :: k, kasave, kovun, kreslt, kwrnsv, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmcosh, fmdiv, fmentr, fmeq, fmeq2, fmexit, &
        fmexp2, fmi2m, fmntr, fmrslt, fmsinh, fmsqr, fmsqrt, fmsub, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      IF (abs(ma(1))>mexpab) THEN
        CALL fmentr('FMTANH',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        namest(ncall) = 'FMTANH'
        IF (ntrace/=0) CALL fmntr(2,ma,ma,1)
        kovun = 0
        IF (ma(1)==mexpov .OR. ma(1)==mexpun) kovun = 1
        ndsave = ndig
        IF (ncall==1) THEN
          k = max(ngrd52-1,2)
          ndig = max(ndig+k,2)
          IF (ndig>ndg2mx) THEN
            kflag = -9
            CALL fmwarn
            ndig = ndsave
            kreslt = 12
            CALL fmrslt(ma,ma,mb,kreslt)
            IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
            ncall = ncall - 1
            RETURN
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 0
        mxsave = mxexp
        mxexp = mxexp2
      END IF

      kwrnsv = kwarn
      kwarn = 0
      ma2 = ma(2)

      macca = ma(0)
      CALL fmeq2(ma,mb,ndsave,ndig,0)
      mb(0) = nint(ndig*alogm2)
      mb(2) = abs(mb(2))
      IF (ma(2)==0) THEN
        CALL fmeq(ma,mb)
        GO TO 10
      END IF

      IF (ma(1)>=1) THEN
        xt = real((ndig+1)/2)*alogmb
        k = int(log(xt)/alogmb)
        IF (ma(1)>k+1) THEN
          CALL fmi2m(1,mb)
          GO TO 10
        ELSE
          x = real(mb(2)*mbase+mb(3)+1)*real(mbase)**int(mb(1)-2)
          IF (x>xt+5.0) THEN
            CALL fmi2m(1,mb)
            GO TO 10
          END IF
        END IF
      END IF
      IF (mb(1)==0 .AND. ndig<50) THEN
        CALL fmexp2(mb,mb)
        CALL fmsqr(mb,mb)
        CALL fmi2m(1,m02)
        CALL fmsub(mb,m02,m03)
        CALL fmadd(mb,m02,m02)
        CALL fmdiv(m03,m02,mb)
        GO TO 10
      END IF
      IF (mb(1)>=0 .AND. mb(2)/=0) THEN
        CALL fmcosh(mb,mb)
        IF (mb(1)>ndig) THEN
          IF (ma2>0) THEN
            CALL fmi2m(1,mb)
            GO TO 10
          ELSE
            CALL fmi2m(-1,mb)
            GO TO 10
          END IF
        END IF
        CALL fmsqr(mb,m03)
        CALL fmi2m(-1,m02)
        CALL fmadd(m03,m02,m03)
        CALL fmsqrt(m03,m03)
        CALL fmdiv(m03,mb,mb)
      ELSE
        CALL fmsinh(mb,mb)
        CALL fmsqr(mb,m03)
        CALL fmi2m(1,m02)
        CALL fmadd(m03,m02,m03)
        CALL fmsqrt(m03,m03)
        CALL fmdiv(mb,m03,mb)
      END IF

!             Round and return.

10    kwarn = kwrnsv
      macmax = nint((ndsave-1)*alogm2+log(real(abs(mb(2))+1))/0.69315)
      mb(0) = min(mb(0),macca,macmax)
      IF (ma2<0 .AND. mb(1)/=munkno) mb(2) = -mb(2)
      CALL fmexit(mb,mb,ndsave,mxsave,kasave,kovun)
      RETURN
    END SUBROUTINE fmtanh
    SUBROUTINE fmtrap(ma)

!  If MA has overflowed or underflowed, replace it by the appropriate
!  symbol.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (ncall<=0) RETURN
      IF (ma(1)>mxexp+1) THEN
        ma(1) = mexpov
        IF (ma(2)>0) THEN
          ma(2) = 1
        ELSE
          ma(2) = -1
        END IF
        ma(0) = nint(ndig*alogm2)
        kflag = -5
      END IF
      IF (ma(1)<-mxexp) THEN
        ma(1) = mexpun
        IF (ma(2)>0) THEN
          ma(2) = 1
        ELSE
          ma(2) = -1
        END IF
        ma(0) = nint(ndig*alogm2)
        kflag = -6
      END IF

      RETURN
    END SUBROUTINE fmtrap
    SUBROUTINE fmulp(ma,mb)

!  MB = The value of one Unit in the Last Place of MA at the current
!       base and precision.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma1
      INTEGER :: j, kwrnsv, n1
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmim, fmmove, fmntr, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      kflag = 0
      ncall = ncall + 1
      namest(ncall) = 'FMULP '
      IF (ntrace/=0) CALL fmntr(2,ma,ma,1)

      ma1 = ma(1)
      n1 = ndig + 1
      DO 10 j = 3, n1
        mwa(j) = 0
10    CONTINUE
      mwa(2) = 1
      IF (ma(2)<0) mwa(2) = -1
      mwa(1) = ma(1) - ndig + 1
      IF (ma(2)==0 .OR. ma(1)>=mexpov) THEN
        CALL fmim(0,mb)
        mb(1) = munkno
        mb(2) = 1
        mb(0) = nint(ndig*alogm2)
        kflag = -4
        IF (ma1/=munkno) CALL fmwarn
      ELSE
        kwrnsv = kwarn
        IF (ma1==mexpun) kwarn = 0
        CALL fmmove(mwa,mb)
        IF (kflag<0) THEN
          namest(ncall) = 'FMULP '
          CALL fmwarn
        END IF
        kwarn = kwrnsv
      END IF
      mb(0) = nint(ndig*alogm2)

      IF (ntrace/=0) CALL fmntr(1,mb,mb,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE fmulp
    SUBROUTINE fmunpk(mp,ma)

!  MP is unpacked and the value returned in MA.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dint, mod
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mp(0:lpack)
! ..
! .. Local Scalars ..
      INTEGER :: j, kp
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: mbase
      INTEGER :: jform1, jform2, kdebug, keswch, kflag, krad, kround, kswide, &
        kw, kwarn, lvltrc, ndig, ntrace
! ..
! .. Common Blocks ..
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      kp = 2
      ma(0) = mp(0)
      ma(1) = mp(1)
      ma(2) = dint(abs(mp(2))/mbase)
      ma(3) = abs(mp(2)) - ma(2)*mbase
      IF (mp(2)<0) ma(2) = -ma(2)
      IF (ndig>=4) THEN
        DO 10 j = 4, ndig, 2
          kp = kp + 1
          ma(j) = dint(mp(kp)/mbase)
          ma(j+1) = mp(kp) - ma(j)*mbase
10      CONTINUE
      END IF
      IF (mod(ndig,2)==1) ma(ndig+1) = dint(mp(kp+1)/mbase)
      RETURN
    END SUBROUTINE fmunpk
    SUBROUTINE fmwarn

!  Called by one of the FM routines to print a warning message
!  if any error condition arises in that routine.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC int
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Local Scalars ..
      INTEGER :: ncs
      CHARACTER (6) :: name
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (kflag>=0 .OR. ncall/=1 .OR. kwarn<=0) RETURN
      ncs = ncall
      name = namest(ncall)
      WRITE (kw,90000) kflag, name

10    ncall = ncall - 1
      IF (ncall>0) THEN
        name = namest(ncall)
        WRITE (kw,90010) name
        GO TO 10
      END IF

      IF (kflag==-1) THEN
        WRITE (kw,90020) ndigmx
      ELSE IF (kflag==-2) THEN
        WRITE (kw,90030) int(mxbase)
      ELSE IF (kflag==-3) THEN
        WRITE (kw,90040)
        WRITE (kw,90050)
      ELSE IF (kflag==-4 .OR. kflag==-7) THEN
        WRITE (kw,90060)
        WRITE (kw,90050)
      ELSE IF (kflag==-5) THEN
        WRITE (kw,90070)
      ELSE IF (kflag==-6) THEN
        WRITE (kw,90080)
      ELSE IF (kflag==-8 .AND. name=='FMOUT ') THEN
        WRITE (kw,90090)
      ELSE IF (kflag==-8 .AND. name=='FMREAD') THEN
        WRITE (kw,90100)
      ELSE IF (kflag==-9) THEN
        WRITE (kw,90110)
        WRITE (kw,90120) ndig, ndg2mx
        WRITE (kw,90050)
      ELSE IF (kflag==-10) THEN
        IF (namest(ncs)=='FMM2SP') THEN
          WRITE (kw,90130)
        ELSE
          WRITE (kw,90140)
        END IF
        WRITE (kw,90150)
      END IF

      ncall = ncs
      IF (kwarn>=2) THEN
        STOP
      END IF
      RETURN
90000 FORMAT (/' Error of type KFLAG =',I3,' in FM package in routine ',A6/)
90010 FORMAT (' called from ',A6)
90020 FORMAT (' NDIG must be between 2 and',I10/)
90030 FORMAT (' MBASE must be between 2 and',I10/)
90040 FORMAT (' An input argument is not a valid FM number.', &
        '  Its exponent is out of range.'/)
90050 FORMAT (' UNKNOWN has been returned.'/)
90060 FORMAT (' Invalid input argument for this routine.'/)
90070 FORMAT (' The result has overflowed.'/)
90080 FORMAT (' The result has underflowed.'/)
90090 FORMAT (' The result array is not big enough to hold the', &
        ' output character string'/' in the current format.'/ &
        ' The result ''***...***'' has been returned.'/)
90100 FORMAT (' The CMBUFF array is not big enough to hold the', &
        ' input character string'/' UNKNOWN has been returned.'/)
90110 FORMAT (' Precision could not be raised enough to provide all', &
        ' requested guard digits.'/)
90120 FORMAT (I23,' digits were requested (NDIG).'/ &
        ' Maximum number of digits currently available',' (NDG2MX) is',I7, &
        '.'/)
90130 FORMAT (' An FM number was too small in magnitude to ', &
        'convert to single precision.'/)
90140 FORMAT (' An FM number was too small in magnitude to ', &
        'convert to double precision.'/)
90150 FORMAT (' Zero has been returned.'/)
    END SUBROUTINE fmwarn
    SUBROUTINE fmwrit(kwrite,ma)

!  Write MA on unit KWRITE.  Multi-line numbers will have '&' as the
!  last nonblank character on all but the last line.  These numbers can
!  then be read easily using FMREAD.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC int, log10, max, min, mod, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kwrite
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: j, jf1sav, jf2sav, k, ksave, l, last, lb, nd, ndsave, nexp
! ..
! .. External Subroutines ..
      EXTERNAL fmeq2, fmout
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      namest(ncall) = 'FMWRIT'
      ndsave = ndig
      ndig = min(ndg2mx,max(ndig+ngrd52,2))

      CALL fmeq2(ma,m01,ndsave,ndig,0)
      ksave = kflag
      nd = int(real(ndig)*log10(real(mbase))) + 1
      IF (nd<2) nd = 2
      nexp = int(2.0*log10(real(mxbase))) + 6
      lb = min(nd+nexp,lmbuff)

      jf1sav = jform1
      jf2sav = jform2
      jform1 = 1
      jform2 = nd + 6

      CALL fmout(m01,cmbuff,lb)

      kflag = ksave
      ndig = ndsave
      jform1 = jf1sav
      jform2 = jf2sav
      last = lb + 1
      DO 10 j = 1, lb
        IF (cmbuff(last-j)/=' ' .OR. j==lb) THEN
          l = last - j
          IF (mod(l,73)/=0) THEN
            WRITE (kwrite,90000) (cmbuff(k),k=1,l)
          ELSE
            IF (l>73) WRITE (kwrite,90000) (cmbuff(k),k=1,l-73)
            WRITE (kwrite,90010) (cmbuff(k),k=l-72,l)
          END IF
          ncall = ncall - 1
          RETURN
        END IF
10    CONTINUE
      ncall = ncall - 1
      RETURN
90000 FORMAT (4X,73A1,' &')
90010 FORMAT (4X,73A1)
    END SUBROUTINE fmwrit

!  Here are the routines that work with packed FM numbers.  All names
!  are the same as unpacked versions with 'FM' replaced by 'FP'.

!  To convert a program using the FM package from unpacked calls to
!  packed calls make these changes to the program:
!  '(0:LUNPCK)' to '(0:LPACK)' in dimensions.
!  'CALL FM' to 'CALL FP'
!  'FMCOMP' to 'FPCOMP'.

    SUBROUTINE fpabs(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmabs, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmabs(mx,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fpabs
    SUBROUTINE fpacos(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmacos, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmacos(mx,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fpacos
    SUBROUTINE fpadd(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmunpk(mb,my)
      CALL fmadd(mx,my,mx)
      CALL fmpack(mx,mc)
      RETURN
    END SUBROUTINE fpadd
    SUBROUTINE fpaddi(ma,l)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: l
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmaddi, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmaddi(mx,l)
      CALL fmpack(mx,ma)
      RETURN
    END SUBROUTINE fpaddi
    SUBROUTINE fpasin(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmasin, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmasin(mx,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fpasin
    SUBROUTINE fpatan(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmatan, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmatan(mx,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fpatan
    SUBROUTINE fpatn2(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmatn2, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmunpk(mb,my)
      CALL fmatn2(mx,my,mx)
      CALL fmpack(mx,mc)
      RETURN
    END SUBROUTINE fpatn2
    SUBROUTINE fpbig(ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmbig, fmpack
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmbig(mx)
      CALL fmpack(mx,ma)
      RETURN
    END SUBROUTINE fpbig
    SUBROUTINE fpchsh(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmchsh, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmchsh(mx,my,mx)
      CALL fmpack(my,mb)
      CALL fmpack(mx,mc)
      RETURN
    END SUBROUTINE fpchsh
    FUNCTION fpcomp(ma,lrel,mb)
      IMPLICIT NONE
! .. Function Return Value ..
      LOGICAL :: fpcomp
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (2) :: lrel
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmunpk(mb,my)
      fpcomp = fmcomp(mx,lrel,my)
      RETURN
    END FUNCTION fpcomp
    SUBROUTINE fpcos(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmcos, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmcos(mx,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fpcos
    SUBROUTINE fpcosh(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmcosh, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmcosh(mx,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fpcosh
    SUBROUTINE fpcssn(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmcssn, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmcssn(mx,my,mx)
      CALL fmpack(my,mb)
      CALL fmpack(mx,mc)
      RETURN
    END SUBROUTINE fpcssn
    SUBROUTINE fpdig(nstack,kst)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kst
! ..
! .. Array Arguments ..
      INTEGER :: nstack(19)
! ..
! .. External Subroutines ..
      EXTERNAL fmdig
! ..

      CALL fmdig(nstack,kst)
      RETURN
    END SUBROUTINE fpdig
    SUBROUTINE fpdim(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmdim, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmunpk(mb,my)
      CALL fmdim(mx,my,mx)
      CALL fmpack(mx,mc)
      RETURN
    END SUBROUTINE fpdim
    SUBROUTINE fpdiv(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmdiv, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmunpk(mb,my)
      CALL fmdiv(mx,my,mx)
      CALL fmpack(mx,mc)
      RETURN
    END SUBROUTINE fpdiv
    SUBROUTINE fpdivi(ma,ival,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmdivi, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmdivi(mx,ival,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fpdivi
    SUBROUTINE fpdp2m(x,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL (KIND(0.0D0)) :: x
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmdp2m, fmpack
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmdp2m(x,mx)
      CALL fmpack(mx,ma)
      RETURN
    END SUBROUTINE fpdp2m
    SUBROUTINE fpdpm(x,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL (KIND(0.0D0)) :: x
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmdpm, fmpack
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmdpm(x,mx)
      CALL fmpack(mx,ma)
      RETURN
    END SUBROUTINE fpdpm
    SUBROUTINE fpeq(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmeq, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmeq(mx,my)
      CALL fmpack(my,mb)
      RETURN
    END SUBROUTINE fpeq
    SUBROUTINE fpequ(ma,mb,nda,ndb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: nda, ndb
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. Local Scalars ..
      INTEGER :: ndasav, ndbsav, ndgsav
! ..
! .. External Subroutines ..
      EXTERNAL fmeq2, fmpack, fmunpk
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: mbase
      INTEGER :: jform1, jform2, kdebug, keswch, kflag, krad, kround, kswide, &
        kw, kwarn, lvltrc, ndig, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..

      ndgsav = ndig
      ndasav = nda
      ndbsav = ndb
      ndig = ndasav
      CALL fmunpk(ma,mx)
      CALL fmeq2(mx,mx,ndasav,ndbsav,1)
      ndig = ndbsav
      CALL fmpack(mx,mb)
      nda = ndasav
      ndb = ndbsav
      ndig = ndgsav
      RETURN
    END SUBROUTINE fpequ
    SUBROUTINE fpexp(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmexp, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmexp(mx,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fpexp
    SUBROUTINE fpform(form,ma,string)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (*) :: form, string
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmform, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmform(form,mx,string)
      RETURN
    END SUBROUTINE fpform
    SUBROUTINE fpfprt(form,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (*) :: form
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmfprt, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmfprt(form,mx)
      RETURN
    END SUBROUTINE fpfprt
    SUBROUTINE fpi2m(ival,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmi2m, fmpack
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmi2m(ival,mx)
      CALL fmpack(mx,ma)
      RETURN
    END SUBROUTINE fpi2m
    SUBROUTINE fpinp(line,ma,la,lb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: la, lb
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
      CHARACTER (1) :: line(lb)
! ..
! .. External Subroutines ..
      EXTERNAL fminp, fmpack
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fminp(line,mx,la,lb)
      CALL fmpack(mx,ma)
      RETURN
    END SUBROUTINE fpinp
    SUBROUTINE fpint(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmint, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmint(mx,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fpint
    SUBROUTINE fpipwr(ma,ival,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmipwr, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmipwr(mx,ival,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fpipwr
    SUBROUTINE fplg10(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmlg10, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmlg10(mx,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fplg10
    SUBROUTINE fpln(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmln, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmln(mx,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fpln
    SUBROUTINE fplni(ival,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmlni, fmpack
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmlni(ival,mx)
      CALL fmpack(mx,ma)
      RETURN
    END SUBROUTINE fplni
    SUBROUTINE fpm2dp(ma,x)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL (KIND(0.0D0)) :: x
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmm2dp, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmm2dp(mx,x)
      RETURN
    END SUBROUTINE fpm2dp
    SUBROUTINE fpm2i(ma,ival)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmm2i, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmm2i(mx,ival)
      RETURN
    END SUBROUTINE fpm2i
    SUBROUTINE fpm2sp(ma,x)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL :: x
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmm2sp, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmm2sp(mx,x)
      RETURN
    END SUBROUTINE fpm2sp
    SUBROUTINE fpmax(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmmax, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmunpk(mb,my)
      CALL fmmax(mx,my,mx)
      CALL fmpack(mx,mc)
      RETURN
    END SUBROUTINE fpmax
    SUBROUTINE fpmin(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmmin, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmunpk(mb,my)
      CALL fmmin(mx,my,mx)
      CALL fmpack(mx,mc)
      RETURN
    END SUBROUTINE fpmin
    SUBROUTINE fpmod(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmmod, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmunpk(mb,my)
      CALL fmmod(mx,my,mx)
      CALL fmpack(mx,mc)
      RETURN
    END SUBROUTINE fpmod
    SUBROUTINE fpmpy(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmmpy, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmunpk(mb,my)
      CALL fmmpy(mx,my,mx)
      CALL fmpack(mx,mc)
      RETURN
    END SUBROUTINE fpmpy
    SUBROUTINE fpmpyi(ma,ival,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmmpyi, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmmpyi(mx,ival,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fpmpyi
    SUBROUTINE fpnint(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmnint, fmpack, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmnint(mx,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fpnint
    SUBROUTINE fpout(ma,line,lb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: lb
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
      CHARACTER (1) :: line(lb)
! ..
! .. External Subroutines ..
      EXTERNAL fmout, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmout(mx,line,lb)
      RETURN
    END SUBROUTINE fpout
    SUBROUTINE fppi(ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack, fmpi
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmpi(mx)
      CALL fmpack(mx,ma)
      RETURN
    END SUBROUTINE fppi
    SUBROUTINE fpprnt(ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmprnt, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmprnt(mx)
      RETURN
    END SUBROUTINE fpprnt
    SUBROUTINE fppwr(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack, fmpwr, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmunpk(mb,my)
      CALL fmpwr(mx,my,mx)
      CALL fmpack(mx,mc)
      RETURN
    END SUBROUTINE fppwr
    SUBROUTINE fpread(kread,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kread
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack, fmread
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmread(kread,mx)
      CALL fmpack(mx,ma)
      RETURN
    END SUBROUTINE fpread
    SUBROUTINE fprpwr(ma,kval,jval,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: jval, kval
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack, fmrpwr, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmrpwr(mx,kval,jval,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fprpwr
    SUBROUTINE fpset(nprec)
! .. Scalar Arguments ..
      INTEGER :: nprec
! ..
! .. External Subroutines ..
      EXTERNAL fmset
! ..

      CALL fmset(nprec)
      RETURN
    END SUBROUTINE fpset
    SUBROUTINE fpsign(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack, fmsign, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmunpk(mb,my)
      CALL fmsign(mx,my,mx)
      CALL fmpack(mx,mc)
      RETURN
    END SUBROUTINE fpsign
    SUBROUTINE fpsin(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack, fmsin, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmsin(mx,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fpsin
    SUBROUTINE fpsinh(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack, fmsinh, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmsinh(mx,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fpsinh
    SUBROUTINE fpsp2m(x,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL :: x
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack, fmsp2m
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmsp2m(x,mx)
      CALL fmpack(mx,ma)
      RETURN
    END SUBROUTINE fpsp2m
    SUBROUTINE fpsqr(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack, fmsqr, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmsqr(mx,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fpsqr
    SUBROUTINE fpsqrt(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack, fmsqrt, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmsqrt(mx,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fpsqrt
    SUBROUTINE fpst2m(string,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (*) :: string
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack, fmst2m
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmst2m(string,mx)
      CALL fmpack(mx,ma)
      RETURN
    END SUBROUTINE fpst2m
    SUBROUTINE fpsub(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack, fmsub, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmunpk(mb,my)
      CALL fmsub(mx,my,mx)
      CALL fmpack(mx,mc)
      RETURN
    END SUBROUTINE fpsub
    SUBROUTINE fptan(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack, fmtan, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmtan(mx,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fptan
    SUBROUTINE fptanh(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack, fmtanh, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmtanh(mx,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fptanh
    SUBROUTINE fpulp(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack, fmulp, fmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmulp(mx,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE fpulp
    SUBROUTINE fpwrit(kwrite,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kwrite
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmunpk, fmwrit
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL fmwrit(kwrite,mx)
      RETURN
    END SUBROUTINE fpwrit

!  The IM routines perform integer multiple-precision arithmetic.




    SUBROUTINE imabs(ma,mb)

!  MB = ABS(MA)

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: kwrnsv, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL imargs, imeq, imntr
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMABS ',1,ma,ma)
      ndsave = ndig
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMABS '
        CALL imntr(2,ma,ma,1)
      END IF

      kflag = 0
      kwrnsv = kwarn
      kwarn = 0
      CALL imeq(ma,mb)
      mb(2) = abs(mb(2))
      kwarn = kwrnsv

      IF (ntrace/=0) CALL imntr(1,mb,mb,1)
      ncall = ncall - 1
      ndig = ndsave
      RETURN
    END SUBROUTINE imabs
    SUBROUTINE imadd(ma,mb,mc)

!  MC = MA + MB

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dint, int, min, nint, sign
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: mda, mdab, mdb
      INTEGER :: j, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmadd2, fmwarn, imargs, imntr
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMADD ',2,ma,mb)
      ndsave = ndig
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMADD '
        CALL imntr(2,ma,mb,2)
      END IF
      kflag = 0

      IF (ma(1)<=2) THEN
        IF (mb(1)>2 .OR. ma(1)<0 .OR. mb(1)<0) GO TO 10
        IF (ma(1)<=1) THEN
          mda = ma(2)
        ELSE IF (ma(2)<0) THEN
          mda = ma(2)*mbase - ma(3)
        ELSE
          mda = ma(2)*mbase + ma(3)
        END IF
        IF (mb(1)<=1) THEN
          mdb = mb(2)
        ELSE IF (mb(2)<0) THEN
          mdb = mb(2)*mbase - mb(3)
        ELSE
          mdb = mb(2)*mbase + mb(3)
        END IF
        mdab = mda + mdb
        IF (abs(mdab)<mbase) THEN
          mc(0) = min(ma(0),mb(0))
          mc(1) = 1
          IF (mdab==0) mc(1) = 0
          mc(2) = mdab
          mc(3) = 0
          IF (mda==0 .OR. mdb==0) kflag = 1
          GO TO 40
        ELSE IF (abs(mdab)<mbase*mbase) THEN
          mc(0) = min(ma(0),mb(0))
          mc(1) = 2
          mc(2) = dint(mdab/mbase)
          mc(3) = abs(mdab-mbase*mc(2))
          IF (mda==0 .OR. mdb==0) kflag = 1
          GO TO 40
        END IF
      END IF

!             Check for special cases.

10    IF (ma(1)>ndg2mx .OR. mb(1)>ndg2mx .OR. ma(1)<0 .OR. mb(1)<0) THEN
        IF (ma(1)==munkno .OR. mb(1)==munkno) THEN
          mc(1) = munkno
          mc(2) = 1
          mc(3) = 0
          mc(0) = nint(ndg2mx*alogm2)
          kflag = -4
          GO TO 50
        END IF
        IF (ma(1)==mexpov) THEN
          mda = 1
          IF ((sign(mda,ma(2))==sign(mda,mb(2))) .OR. (mb(2)==0)) THEN
            mc(0) = ma(0)
            mc(1) = ma(1)
            mc(2) = ma(2)
            mc(3) = ma(3)
            kflag = -5
            GO TO 50
          ELSE
            mc(0) = nint(ndg2mx*alogm2)
            mc(1) = munkno
            mc(2) = 1
            mc(3) = 0
            kflag = -4
            namest(ncall) = 'IMADD '
            CALL fmwarn
            GO TO 50
          END IF
        END IF
        IF (mb(1)==mexpov) THEN
          mda = 1
          IF ((sign(mda,mb(2))==sign(mda,ma(2))) .OR. (ma(2)==0)) THEN
            mc(0) = mb(0)
            mc(1) = mb(1)
            mc(2) = mb(2)
            mc(3) = mb(3)
            kflag = -5
            GO TO 50
          ELSE
            mc(0) = nint(ndg2mx*alogm2)
            mc(1) = munkno
            mc(2) = 1
            mc(3) = 0
            kflag = -4
            namest(ncall) = 'IMADD '
            CALL fmwarn
            GO TO 50
          END IF
        END IF
        mc(1) = munkno
        mc(2) = 1
        mc(3) = 0
        mc(0) = nint(ndg2mx*alogm2)
        kflag = -4
        namest(ncall) = 'IMADD '
        CALL fmwarn
        GO TO 50
      END IF

      IF (ma(1)>mb(1)) THEN
        ndig = int(ma(1)) + 1
        IF (ndig<2 .OR. ndig>ndg2mx) ndig = 2
        ma(ndig+1) = 0
        DO 20 j = int(mb(1)) + 2, ndig + 1
          mb(j) = 0
20      CONTINUE
      ELSE
        ndig = int(mb(1)) + 1
        IF (ndig<2 .OR. ndig>ndg2mx) ndig = 2
        mb(ndig+1) = 0
        DO 30 j = int(ma(1)) + 2, ndig + 1
          ma(j) = 0
30      CONTINUE
      END IF

      CALL fmadd2(ma,mb,mc)

40    IF (mc(1)>ndigmx) THEN
        IF (ncall==1 .OR. mc(1)>ndg2mx) THEN
          mc(0) = nint(ndg2mx*alogm2)
          mc(1) = mexpov
          IF (mc(2)>0) THEN
            mc(2) = 1
          ELSE
            mc(2) = -1
          END IF
          mc(3) = 0
          kflag = -5
          namest(ncall) = 'IMADD '
          CALL fmwarn
        END IF
      END IF

50    IF (ntrace/=0) CALL imntr(1,mc,mc,1)
      ncall = ncall - 1
      ndig = ndsave
      RETURN
    END SUBROUTINE imadd
    SUBROUTINE imargs(kroutn,nargs,ma,mb)

!  Check the input arguments to a routine for special cases.

!  KROUTN - Name of the subroutine that was called
!  NARGS  - The number of input arguments (1 or 2)
!  MA     - First input argument
!  MB     - Second input argument (if NARGS is 2)

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, int, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: nargs
      CHARACTER (6) :: kroutn
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: mbs
      INTEGER :: j, kwrnsv, last
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      kflag = -4
      IF (ma(1)==munkno) RETURN
      IF (nargs==2) THEN
        IF (mb(1)==munkno) RETURN
      END IF
      IF (mblogs/=mbase) CALL fmcons
      kflag = 0

!             Check the validity of parameters.

      IF (ncall>1 .AND. kdebug==0) RETURN
      namest(ncall) = kroutn

!             Check MBASE.

      IF (mbase<2 .OR. mbase>mxbase) THEN
        kflag = -2
        CALL fmwarn
        mbs = mbase
        IF (mbase<2) mbase = 2
        IF (mbase>mxbase) mbase = mxbase
        WRITE (kw,90000) int(mbs), int(mbase)
        CALL fmcons
        RETURN
      END IF

!             Check exponent range.

      IF (ma(1)>lunpck .OR. ma(1)<0) THEN
        IF (abs(ma(1))/=mexpov .OR. abs(ma(2))/=1) THEN
          kflag = -3
          CALL fmwarn
          ma(0) = nint(ndg2mx*alogm2)
          ma(1) = munkno
          ma(2) = 1
          ma(3) = 0
          RETURN
        END IF
      END IF
      IF (nargs==2) THEN
        IF (mb(1)>lunpck .OR. mb(1)<0) THEN
          IF (abs(mb(1))/=mexpov .OR. abs(mb(2))/=1) THEN
            kflag = -3
            CALL fmwarn
            mb(0) = nint(ndg2mx*alogm2)
            mb(1) = munkno
            mb(2) = 1
            mb(3) = 0
            RETURN
          END IF
        END IF
      END IF

!             Check for properly normalized digits in the
!             input arguments.

      IF (abs(ma(1)-int(ma(1)))/=0) kflag = 1
      IF (ma(2)<=(-mbase) .OR. ma(2)>=mbase .OR. abs(ma(2)-int(ma(2)))/=0) &
        kflag = 2
      IF (kdebug==0) GO TO 20
      last = int(ma(1)) + 1
      IF (ma(1)>lunpck) last = 3
      DO 10 j = 3, last
        IF (ma(j)<0 .OR. ma(j)>=mbase .OR. abs(ma(j)-int(ma(j)))/=0) THEN
          kflag = j
          GO TO 20
        END IF
10    CONTINUE
20    IF (kflag/=0) THEN
        j = kflag
        mbs = ma(j)
        kflag = -4
        kwrnsv = kwarn
        IF (kwarn>=2) kwarn = 1
        CALL fmwarn
        kwarn = kwrnsv
        IF (kwarn>=1) THEN
          WRITE (kw,*) ' First invalid array element:  MA(', j, ') = ', mbs
        END IF
        ma(0) = nint(ndg2mx*alogm2)
        ma(1) = munkno
        ma(2) = 1
        ma(3) = 0
        IF (kwarn>=2) THEN
          STOP
        END IF
        RETURN
      END IF
      IF (nargs==2) THEN
        IF (abs(mb(1)-int(mb(1)))/=0) kflag = 1
        IF (mb(2)<=(-mbase) .OR. mb(2)>=mbase .OR. abs(mb(2)-int(mb(2)))/=0) &
          kflag = 2
        IF (kdebug==0) GO TO 40
        last = int(mb(1)) + 1
        IF (mb(1)>lunpck) last = 3
        DO 30 j = 3, last
          IF (mb(j)<0 .OR. mb(j)>=mbase .OR. abs(mb(j)-int(mb(j)))/=0) THEN
            kflag = j
            GO TO 40
          END IF
30      CONTINUE
40      IF (kflag/=0) THEN
          j = kflag
          mbs = mb(j)
          kflag = -4
          kwrnsv = kwarn
          IF (kwarn>=2) kwarn = 1
          CALL fmwarn
          kwarn = kwrnsv
          IF (kwarn>=1) THEN
            WRITE (kw,*) ' First invalid array element:  MB(', j, ') = ', mbs
          END IF
          mb(0) = nint(ndg2mx*alogm2)
          mb(1) = munkno
          mb(2) = 1
          mb(3) = 0
          IF (kwarn>=2) THEN
            STOP
          END IF
          RETURN
        END IF
      END IF
      RETURN
90000 FORMAT (' MBASE was',I10,'.  It has been changed to',I10,'.')
    END SUBROUTINE imargs
    SUBROUTINE imbig(ma)

!     MA = The biggest representable IM integer.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: j
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, imntr
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      namest(ncall) = 'IMBIG '

      IF (mblogs/=mbase) CALL fmcons
      kflag = 0
      DO 10 j = 2, ndigmx + 1
        ma(j) = mbase - 1
10    CONTINUE
      ma(1) = ndigmx
      ma(0) = nint(ndigmx*alogm2)

      IF (ntrace/=0) CALL imntr(1,ma,ma,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE imbig
    FUNCTION imcomp(ma,lrel,mb)

!  Logical comparison of FM numbers MA and MB.

!  LREL is a CHARACTER *2 description of the comparison to be done:
!  LREL = 'EQ' returns IMCOMP = .TRUE. if MA.EQ.MB
!       = 'NE', 'GE', 'GT', 'LE', 'LT' also work like a logical IF.

!  Some compilers object to functions with side effects such as
!  changing KFLAG or other common variables.  Blocks of code that
!  modify common are identified by:
!      C                                                 DELETE START
!        ...
!      C                                                 DELETE STOP
!  These may be removed or commented out to produce a function without
!  side effects.  This disables trace printing in IMCOMP, and error
!  codes are not returned in KFLAG.

      IMPLICIT NONE

! .. Function Return Value ..
      LOGICAL :: imcomp
! ..
! .. Intrinsic Functions ..
      INTRINSIC abs, int, max
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (2) :: lrel
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: j, jcomp, ndsave, nlast, ntrsav
      CHARACTER (2) :: jrel
! ..
! .. External Subroutines ..
      EXTERNAL imargs, imntrj, imprnt
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
!                                                 DELETE START
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMCOMP',2,ma,mb)
      namest(ncall) = 'IMCOMP'

      IF (ncall<=lvltrc .AND. abs(ntrace)>=2) THEN
        WRITE (kw,90000)
        ndsave = ndig
        IF (ntrace>0) THEN
          CALL imprnt(ma)
          WRITE (kw,90010) lrel
          CALL imprnt(mb)
        ELSE
          ndig = max(2,int(ma(1)))
          IF (ndig>ndg2mx) ndig = 2
          IF (ma(1)<=1) ma(3) = 0
          ntrsav = ntrace
          IF (ntrace<-2) ntrace = -2
          CALL imntrj(ma,ndig)
          WRITE (kw,90010) lrel
          ndig = max(2,int(mb(1)))
          IF (ndig>ndg2mx) ndig = 2
          IF (mb(1)<=1) mb(3) = 0
          CALL imntrj(mb,ndig)
          ntrace = ntrsav
        END IF
        ndig = ndsave
      END IF
!                                                 DELETE STOP

!             JCOMP will be 1 if MA.GT.MB
!                           2 if MA.EQ.MB
!                           3 if MA.LT.MB

!             Check for special cases.

      jrel = lrel
      IF (lrel/='EQ' .AND. lrel/='NE' .AND. lrel/='LT' .AND. lrel/='GT' .AND. &
          lrel/='LE' .AND. lrel/='GE') THEN
        IF (lrel=='eq') THEN
          jrel = 'EQ'
        ELSE IF (lrel=='ne') THEN
          jrel = 'NE'
        ELSE IF (lrel=='lt') THEN
          jrel = 'LT'
        ELSE IF (lrel=='gt') THEN
          jrel = 'GT'
        ELSE IF (lrel=='le') THEN
          jrel = 'LE'
        ELSE IF (lrel=='ge') THEN
          jrel = 'GE'
        ELSE
          imcomp = .FALSE.
!                                                 DELETE START
          kflag = -4
          IF (ncall/=1 .OR. kwarn<=0) GO TO 30
!                                                 DELETE STOP
          IF (kwarn<=0) GO TO 30
          WRITE (kw,90020) lrel
          IF (kwarn>=2) THEN
            STOP
          END IF
          GO TO 30
        END IF
      END IF

      IF (ma(1)==munkno .OR. mb(1)==munkno) THEN
        imcomp = .FALSE.
!                                                 DELETE START
        kflag = -4
!                                                 DELETE STOP
        GO TO 30
      END IF

      IF (abs(ma(1))==mexpov .AND. ma(1)==mb(1) .AND. ma(2)==mb(2)) THEN
        imcomp = .FALSE.
!                                                 DELETE START
        kflag = -4
        IF (ncall/=1 .OR. kwarn<=0) GO TO 30
!                                                 DELETE STOP
        IF (kwarn<=0) GO TO 30
        WRITE (kw,90030)
        IF (kwarn>=2) THEN
          STOP
        END IF
        GO TO 30
      END IF

!             Check for zero.

!                                                 DELETE START
      kflag = 0
!                                                 DELETE STOP
      IF (ma(2)==0) THEN
        jcomp = 2
        IF (mb(2)<0) jcomp = 1
        IF (mb(2)>0) jcomp = 3
        GO TO 20
      END IF
      IF (mb(2)==0) THEN
        jcomp = 1
        IF (ma(2)<0) jcomp = 3
        GO TO 20
      END IF
!             Check for opposite signs.

      IF (ma(2)>0 .AND. mb(2)<0) THEN
        jcomp = 1
        GO TO 20
      END IF
      IF (mb(2)>0 .AND. ma(2)<0) THEN
        jcomp = 3
        GO TO 20
      END IF

!             See which one is larger in absolute value.

      IF (ma(1)>mb(1)) THEN
        jcomp = 1
        GO TO 20
      END IF
      IF (mb(1)>ma(1)) THEN
        jcomp = 3
        GO TO 20
      END IF
      nlast = int(ma(1)) + 1
      IF (nlast>ndg2mx+1) nlast = 2

      DO 10 j = 2, nlast
        IF (abs(ma(j))>abs(mb(j))) THEN
          jcomp = 1
          GO TO 20
        END IF
        IF (abs(mb(j))>abs(ma(j))) THEN
          jcomp = 3
          GO TO 20
        END IF
10    CONTINUE

      jcomp = 2

!             Now match the JCOMP value to the requested comparison.

20    IF (jcomp==1 .AND. ma(2)<0) THEN
        jcomp = 3
      ELSE IF (jcomp==3 .AND. mb(2)<0) THEN
        jcomp = 1
      END IF

      imcomp = .FALSE.
      IF (jcomp==1 .AND. (jrel=='GT' .OR. jrel=='GE' .OR. jrel=='NE')) &
        imcomp = .TRUE.

      IF (jcomp==2 .AND. (jrel=='EQ' .OR. jrel=='GE' .OR. jrel=='LE')) &
        imcomp = .TRUE.

      IF (jcomp==3 .AND. (jrel=='NE' .OR. jrel=='LT' .OR. jrel=='LE')) &
        imcomp = .TRUE.

30    CONTINUE
!                                                 DELETE START
      IF (ntrace/=0) THEN
        IF (ncall<=lvltrc .AND. abs(ntrace)>=1) THEN
          IF (kflag==0) THEN
            WRITE (kw,90040) ncall, int(mbase)
          ELSE
            WRITE (kw,90050) ncall, int(mbase), kflag
          END IF
          IF (imcomp) THEN
            WRITE (kw,90060)
          ELSE
            WRITE (kw,90070)
          END IF
        END IF
      END IF
      ncall = ncall - 1
!                                                 DELETE STOP
      RETURN
90000 FORMAT (' Input to IMCOMP')
90010 FORMAT (7X,'.',A2,'.')
90020 FORMAT (/' Error of type KFLAG = -4 in FM package in', &
        ' routine IMCOMP'//1X,A,' is not one of the six', &
        ' recognized comparisons.'//' .FALSE. has been',' returned.'/)
90030 FORMAT (/' Error of type KFLAG = -4 in FM package in routine', &
        ' IMCOMP'//' Two numbers in the same overflow', &
        ' category cannot be compared.'//' .FALSE. has been returned.'/)
90040 FORMAT (' IMCOMP',15X,'Call level =',I2,5X,'MBASE =',I10)
90050 FORMAT (' IMCOMP',6X,'Call level =',I2,4X,'MBASE =',I10,4X,'KFLAG =',I3)
90060 FORMAT (7X,'.TRUE.')
90070 FORMAT (7X,'.FALSE.')
    END FUNCTION imcomp
    SUBROUTINE imdim(ma,mb,mc)

!  MC = DIM(MA,MB)

!  Positive difference.  MC = MA - MB  if MA.GE.MB,
!                           = 0        otherwise.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: kovfl
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
      EXTERNAL fmwarn, imargs, imntr, imsub
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      kflag = 0
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMDIM ',2,ma,mb)
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMDIM '
        CALL imntr(2,ma,mb,2)
      END IF

      IF (ma(1)==munkno .OR. mb(1)==munkno) THEN
        mc(1) = munkno
        mc(2) = 1
        mc(3) = 0
        mc(0) = nint(ndg2mx*alogm2)
        kflag = -4
        GO TO 10
      END IF
      IF (ma(1)<0 .OR. mb(1)<0) THEN
        mc(1) = munkno
        mc(2) = 1
        mc(3) = 0
        mc(0) = nint(ndg2mx*alogm2)
        kflag = -4
        namest(ncall) = 'IMDIM '
        CALL fmwarn
        GO TO 10
      END IF
      kovfl = 0
      IF (ma(1)==mexpov .OR. mb(1)==mexpov) THEN
        kovfl = 1
        IF (ma(1)==mexpov .AND. mb(1)==mexpov .AND. ma(2)==mb(2)) THEN
          mc(1) = munkno
          mc(2) = 1
          mc(3) = 0
          mc(0) = nint(ndg2mx*alogm2)
          kflag = -4
          namest(ncall) = 'IMDIM '
          CALL fmwarn
          GO TO 10
        END IF
      END IF

      IF (imcomp(ma,'GE',mb)) THEN
        CALL imsub(ma,mb,mc)
        IF (kflag==1) kflag = 0
      ELSE
        mc(1) = 0
        mc(2) = 0
        mc(3) = 0
        mc(0) = nint(ndg2mx*alogm2)
      END IF

      IF (mc(1)>ndigmx) THEN
        IF (mc(1)==munkno) THEN
          kflag = -4
          namest(ncall) = 'IMDIM '
          CALL fmwarn
        ELSE IF (ncall==1 .OR. mc(1)>ndg2mx) THEN
          mc(0) = nint(ndg2mx*alogm2)
          mc(1) = mexpov
          IF (mc(2)>0) THEN
            mc(2) = 1
          ELSE
            mc(2) = -1
          END IF
          mc(3) = 0
          kflag = -5
          namest(ncall) = 'IMDIM '
          IF (kovfl/=1) CALL fmwarn
        END IF
      END IF

10    IF (ntrace/=0) CALL imntr(1,mc,mc,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE imdim
    SUBROUTINE imdiv(ma,mb,mc)

!  MC = INT(MA/MB)

!  Use IMDIVR if both INT(MA/MB) and MOD(MA,MB) are needed.

      IMPLICIT NONE

!             Scratch array usage during IMDIV:   M01 - M03

! .. Intrinsic Functions ..
      INTRINSIC nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmwarn, imargs, imdivr, imntr
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMDIV ',2,ma,mb)
      kflag = 0
      ndsave = ndig
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMDIV '
        CALL imntr(2,ma,mb,2)
      END IF

      IF (ma(1)==munkno .OR. mb(1)==munkno) THEN
        mc(1) = munkno
        mc(2) = 1
        mc(3) = 0
        mc(0) = nint(ndg2mx*alogm2)
        kflag = -4
        GO TO 10
      END IF

      CALL imdivr(ma,mb,mc,m03)

      IF (mc(1)==munkno) THEN
        kflag = -4
        namest(ncall) = 'IMDIV '
        CALL fmwarn
      END IF

10    IF (ntrace/=0) CALL imntr(1,mc,mc,1)
      ncall = ncall - 1
      ndig = ndsave
      RETURN
    END SUBROUTINE imdiv
    SUBROUTINE imdivi(ma,idiv,mb)

!  MB = INT(MA/IDIV)

!  Use IMDVIR if both INT(MA/IDIV) and MOD(MA,IDIV) are needed.

      IMPLICIT NONE

!             Scratch array usage during IMDIVI:   M01 - M03

! .. Intrinsic Functions ..
      INTRINSIC nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: idiv
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: irem, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmwarn, imargs, imdvir, imntr, imntri
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMDIVI',1,ma,ma)
      kflag = 0
      ndsave = ndig
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMDIVI'
        CALL imntr(2,ma,ma,1)
        CALL imntri(2,idiv,0)
      END IF

      IF (ma(1)==munkno) THEN
        mb(1) = munkno
        mb(2) = 1
        mb(3) = 0
        mb(0) = nint(ndg2mx*alogm2)
        kflag = -4
        GO TO 10
      END IF

      CALL imdvir(ma,idiv,mb,irem)

      IF (mb(1)==munkno) THEN
        kflag = -4
        namest(ncall) = 'IMDIVI'
        CALL fmwarn
      END IF

10    IF (ntrace/=0) CALL imntr(1,mb,mb,1)
      ncall = ncall - 1
      ndig = ndsave
      RETURN
    END SUBROUTINE imdivi
    SUBROUTINE imdivr(ma,mb,mc,md)

!  MC = INT(MA / MB),    MD = Remainder from the division.

      IMPLICIT NONE

!             Scratch array usage during IMDIVR:   M01 - M02

! .. Intrinsic Functions ..
      INTRINSIC abs, dble, dint, int, max, min, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, ma2p, macca, maccb, maxmwa, mb1, mb2, mb2p, mbm1, &
        mcarry, mda, mdab, mdb, mdr, mkt, mlmax, mqd
      REAL (KIND(0.0D0)) :: xb, xbase, xbr, xmwa
      INTEGER :: j, jb, jl, k, ka, kb, kl, kltflg, kptmwa, lcrrct, na1, nb1, &
        ndsave, nguard, nl, nmbwds, ntrsav
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmwarn, imadd, imargs, imeq, imi2m, imntr, imntrj, &
        imprnt, imsub
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMDIVR',2,ma,mb)
      ndsave = ndig
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMDIVR'
        CALL imntr(2,ma,mb,2)
      END IF
      kflag = 0
      ntrsav = ntrace
      ntrace = 0
      IF (mblogs/=mbase) CALL fmcons

!             Check for special cases.

      IF (mb(1)==1 .AND. ma(1)/=munkno) THEN
        IF (mb(2)==1) THEN
          CALL imeq(ma,mc)
          md(1) = 0
          md(2) = 0
          md(3) = 0
          md(0) = nint(ndg2mx*alogm2)
          GO TO 260
        ELSE IF (mb(2)==-1) THEN
          CALL imeq(ma,mc)
          IF (mc(1)/=munkno) mc(2) = -mc(2)
          md(1) = 0
          md(2) = 0
          md(3) = 0
          md(0) = nint(ndg2mx*alogm2)
          GO TO 260
        END IF
      END IF
      IF (ma(1)<mb(1) .AND. mb(1)/=munkno) GO TO 10
      IF (ma(1)>ndg2mx .OR. mb(1)>ndg2mx .OR. ma(1)<0 .OR. mb(1)<0 .OR. &
          mb(2)==0) THEN
        kflag = -4
        IF (ma(1)/=munkno .AND. mb(1)/=munkno) THEN
          namest(ncall) = 'IMDIVR'
          CALL fmwarn
        END IF
        mc(1) = munkno
        mc(2) = 1
        mc(3) = 0
        mc(0) = nint(ndg2mx*alogm2)
        md(1) = munkno
        md(2) = 1
        md(3) = 0
        md(0) = nint(ndg2mx*alogm2)
        GO TO 260
      END IF
      IF (ma(1)<=2) THEN
        IF (mb(1)>2) GO TO 10
        IF (mb(2)==0) GO TO 10
        IF (ma(1)<=1) THEN
          mda = ma(2)
        ELSE IF (ma(2)<0) THEN
          mda = ma(2)*mbase - ma(3)
        ELSE
          mda = ma(2)*mbase + ma(3)
        END IF
        IF (mb(1)<=1) THEN
          mdb = mb(2)
        ELSE IF (mb(2)<0) THEN
          mdb = mb(2)*mbase - mb(3)
        ELSE
          mdb = mb(2)*mbase + mb(3)
        END IF
        mdab = dint(mda/mdb)
        mdr = mda - mdab*mdb
        IF (abs(mdab)<mbase) THEN
          mc(0) = min(ma(0),mb(0))
          mc(1) = 1
          IF (mdab==0) mc(1) = 0
          mc(2) = mdab
          mc(3) = 0
        ELSE IF (abs(mdab)<mbase*mbase) THEN
          mc(0) = min(ma(0),mb(0))
          mc(1) = 2
          mc(2) = dint(mdab/mbase)
          mc(3) = abs(mdab-mbase*mc(2))
        ELSE
          GO TO 10
        END IF
        IF (abs(mdr)<mbase) THEN
          md(0) = mc(0)
          md(1) = 1
          IF (mdr==0) md(1) = 0
          md(2) = mdr
          md(3) = 0
          GO TO 260
        ELSE IF (abs(mdr)<mbase*mbase) THEN
          md(0) = mc(0)
          md(1) = 2
          md(2) = dint(mdr/mbase)
          md(3) = abs(mdr-mbase*md(2))
          GO TO 260
        END IF
      END IF

10    kltflg = 0
      ma2 = ma(2)
      mb2 = mb(2)
      kl = int(mb(1))
      IF (kl>ndg2mx) kl = 2
      DO 20 j = 0, kl + 1
        m01(j) = mb(j)
20    CONTINUE
      m01(2) = abs(m01(2))
      IF (kl==1) m01(3) = 0
      IF (ma(1)==m01(1) .AND. abs(ma(2))<=m01(2)) THEN
        ma(2) = abs(ma(2))
        IF (imcomp(ma,'EQ',m01)) THEN
          kltflg = 2
        ELSE IF (imcomp(ma,'LT',m01)) THEN
          kltflg = 1
        END IF
        ma(2) = ma2
      END IF
      IF (ma(1)<mb(1) .OR. kltflg>=1) THEN
        IF (kltflg/=2) THEN
          CALL imeq(ma,md)
          md(2) = abs(md(2))
          CALL imi2m(0,mc)
        ELSE
          CALL imi2m(1,mc)
          CALL imi2m(0,md)
        END IF
        GO TO 250
      END IF

      ndig = int(ma(1))
      IF (ndig<2) ndig = 2

      macca = ma(0)
      maccb = mb(0)

!             NGUARD is the number of guard digits used.

      nguard = 1
      ma2p = abs(ma(2))
      mb2p = abs(mb(2))
      na1 = int(ma(1)) + 1
      nb1 = int(mb(1)) + 1

!             Copy MA into the working array.

      DO 30 j = 3, na1
        mwa(j+1) = ma(j)
30    CONTINUE
      mwa(1) = ma(1) - mb(1) + 1
      mwa(2) = 0
      nl = na1 + nguard + 3
      DO 40 j = na1 + 2, nl
        mwa(j) = 0
40    CONTINUE

!             Save the sign of MA and MB and then work only with
!             positive numbers.

      ma2 = ma(2)
      mb1 = mb(1)
      mb2 = mb(2)
      ma(2) = ma2p
      mwa(3) = ma(2)
      mb(1) = 0
      mb(2) = mb2p

!             NMBWDS is the number of words of MB used to
!             compute the estimated quotient digit MQD.

      nmbwds = 4
      IF (mbase<100) nmbwds = 7

!             XB is an approximation of MB used in
!             estimating the quotient digits.

      xbase = dble(mbase)
      xb = 0
      jl = nmbwds
      IF (jl<=nb1) THEN
        DO 50 j = 2, jl
          xb = xb*xbase + dble(mb(j))
50      CONTINUE
      ELSE
        DO 60 j = 2, jl
          IF (j<=nb1) THEN
            xb = xb*xbase + dble(mb(j))
          ELSE
            xb = xb*xbase
          END IF
60      CONTINUE
      END IF
      IF (jl+1<=nb1) xb = xb + dble(mb(jl+1))/xbase
      xbr = 1.0D0/xb

!             MLMAX determines when to normalize all of MWA.

      mbm1 = mbase - 1
      mlmax = maxint/mbm1
      mkt = intmax - mbase
      mlmax = min(mlmax,mkt)

!             MAXMWA is an upper bound on the size of values in MWA
!             divided by MBASE-1.  It is used to determine whether
!             normalization can be postponed.

      maxmwa = 0

!             KPTMWA points to the next digit in the quotient.

      kptmwa = 2

!             This is the start of the division loop.

!             XMWA is an approximation of the active part of MWA
!             used in estimating quotient digits.

70    kl = kptmwa + nmbwds - 1
      IF (kl<=nl) THEN
        xmwa = ((dble(mwa(kptmwa))*xbase+dble(mwa(kptmwa+1)))*xbase+dble(mwa( &
          kptmwa+2)))*xbase + dble(mwa(kptmwa+3))
        DO 80 j = kptmwa + 4, kl
          xmwa = xmwa*xbase + dble(mwa(j))
80      CONTINUE
      ELSE
        xmwa = dble(mwa(kptmwa))
        DO 90 j = kptmwa + 1, kl
          IF (j<=nl) THEN
            xmwa = xmwa*xbase + dble(mwa(j))
          ELSE
            xmwa = xmwa*xbase
          END IF
90      CONTINUE
      END IF

!             MQD is the estimated quotient digit.

      mqd = dint(xmwa*xbr)
      IF (mqd<0) mqd = mqd - 1

      IF (mqd>0) THEN
        maxmwa = maxmwa + mqd
      ELSE
        maxmwa = maxmwa - mqd
      END IF

!             See if MWA must be normalized.

      ka = kptmwa + 1
      kb = ka + int(mb1) - 1
      IF (maxmwa>=mlmax) THEN
        DO 100 j = kb, ka, -1
          IF (mwa(j)<0) THEN
            mcarry = int((-mwa(j)-1)/mbase) + 1
            mwa(j) = mwa(j) + mcarry*mbase
            mwa(j-1) = mwa(j-1) - mcarry
          ELSE IF (mwa(j)>=mbase) THEN
            mcarry = -int(mwa(j)/mbase)
            mwa(j) = mwa(j) + mcarry*mbase
            mwa(j-1) = mwa(j-1) - mcarry
          END IF
100     CONTINUE
        xmwa = 0
        IF (kl<=nl) THEN
          DO 110 j = kptmwa, kl
            xmwa = xmwa*xbase + dble(mwa(j))
110       CONTINUE
        ELSE
          DO 120 j = kptmwa, kl
            IF (j<=nl) THEN
              xmwa = xmwa*xbase + dble(mwa(j))
            ELSE
              xmwa = xmwa*xbase
            END IF
120       CONTINUE
        END IF
        mqd = dint(xmwa*xbr)
        IF (mqd<0) mqd = mqd - 1
        IF (mqd>0) THEN
          maxmwa = mqd
        ELSE
          maxmwa = -mqd
        END IF
      END IF

!             Subtract MQD*MB from MWA.

      jb = ka - 2
      IF (mqd/=0) THEN

!             Major (Inner Loop)

        DO 130 j = ka, kb
          mwa(j) = mwa(j) - mqd*mb(j-jb)
130     CONTINUE
      END IF

      mwa(ka) = mwa(ka) + mwa(ka-1)*mbase
      mwa(kptmwa) = mqd

      kptmwa = kptmwa + 1
      IF (kptmwa-2<mwa(1)) GO TO 70

!             Final normalization.

      kptmwa = kptmwa - 1
      DO 140 j = kptmwa, 3, -1
        IF (mwa(j)<0) THEN
          mcarry = int((-mwa(j)-1)/mbase) + 1
          mwa(j) = mwa(j) + mcarry*mbase
          mwa(j-1) = mwa(j-1) - mcarry
        ELSE IF (mwa(j)>=mbase) THEN
          mcarry = -int(mwa(j)/mbase)
          mwa(j) = mwa(j) + mcarry*mbase
          mwa(j-1) = mwa(j-1) - mcarry
        END IF
140   CONTINUE

      lcrrct = 0
150   DO 160 j = kptmwa + int(mb1), kptmwa + 2, -1
        IF (mwa(j)<0) THEN
          mcarry = int((-mwa(j)-1)/mbase) + 1
          mwa(j) = mwa(j) + mcarry*mbase
          mwa(j-1) = mwa(j-1) - mcarry
        ELSE IF (mwa(j)>=mbase) THEN
          mcarry = -int(mwa(j)/mbase)
          mwa(j) = mwa(j) + mcarry*mbase
          mwa(j-1) = mwa(j-1) - mcarry
        END IF
160   CONTINUE

!             Due to rounding, the remainder may not be between
!             0 and ABS(MB) here.  Correct if necessary.

      IF (mwa(ka)<0) THEN
        lcrrct = lcrrct - 1
        DO 170 j = ka, kb
          mwa(j) = mwa(j) + mb(j-jb)
170     CONTINUE
        GO TO 150
      ELSE IF (mwa(ka)>=mbase) THEN
        lcrrct = lcrrct + 1
        DO 180 j = ka, kb
          mwa(j) = mwa(j) - mb(j-jb)
180     CONTINUE
        GO TO 150
      END IF

      ma(2) = ma2
      mb(1) = mb1
      mb(2) = mb2

      IF (mwa(2)/=0 .OR. kptmwa==2) THEN
        DO 190 j = 1, int(mwa(1)) + 1
          mc(j) = mwa(j)
190     CONTINUE
      ELSE
        DO 200 j = 3, int(mwa(1)) + 1
          mc(j-1) = mwa(j)
200     CONTINUE
        IF (mc(2)/=0) THEN
          mc(1) = mwa(1) - 1
        ELSE
          mc(1) = 0
        END IF
      END IF
      IF (mc(1)<=1) mc(3) = 0
      mc(0) = min(macca,maccb)

      IF (mwa(kptmwa+1)/=0) THEN
        DO 210 j = 1, int(mb1)
          md(j+1) = mwa(kptmwa+j)
210     CONTINUE
        md(1) = mb1
      ELSE
        DO 230 j = 1, int(mb1)
          IF (mwa(kptmwa+j)/=0) THEN
            DO 220 k = j, int(mb1)
              md(k-j+2) = mwa(kptmwa+k)
220         CONTINUE
            md(1) = mb1 + 1 - j
            GO TO 240
          END IF
230     CONTINUE
        md(1) = 0
        md(2) = 0
      END IF
240   IF (md(1)<=1) md(3) = 0
      md(0) = min(macca,maccb)

!             If the remainder had to be corrected, make the
!             corresponding adjustment in the quotient.

      IF (md(1)>m01(1) .OR. (md(1)==m01(1) .AND. abs(md(2))>=m01(2))) THEN
        IF (imcomp(md,'GE',m01)) THEN
          CALL imsub(md,m01,md)
          lcrrct = lcrrct + 1
        END IF
      END IF
      IF (lcrrct/=0) THEN
        CALL imi2m(lcrrct,m02)
        CALL imadd(m02,mc,mc)
      END IF

250   IF (ma2<0 .AND. mb2>0) THEN
        IF (mc(1)/=munkno) mc(2) = -mc(2)
        IF (md(1)/=munkno) md(2) = -md(2)
      ELSE IF (ma2>0 .AND. mb2<0) THEN
        IF (mc(1)/=munkno) mc(2) = -mc(2)
      ELSE IF (ma2<0 .AND. mb2<0) THEN
        IF (md(1)/=munkno) md(2) = -md(2)
      END IF

260   ntrace = ntrsav
      IF (ntrace/=0) THEN
        CALL imntr(1,mc,mc,1)
        IF (abs(ntrace)>=1 .AND. ncall<=lvltrc) THEN
          IF (ntrace<0) THEN
            ndig = max(2,int(md(1)))
            IF (ndig>ndg2mx) ndig = 2
            IF (md(1)<=1) md(3) = 0
            ntrsav = ntrace
            IF (ntrace<-2) ntrace = -2
            CALL imntrj(md,ndig)
            ntrace = ntrsav
          ELSE
            CALL imprnt(md)
          END IF
        END IF
      END IF
      ncall = ncall - 1
      ndig = ndsave
      RETURN
    END SUBROUTINE imdivr
    SUBROUTINE imdvir(ma,idiv,mb,irem)

!  MB = INT(MA / IDIV),    IREM = Remainder from the division.

!  Division by a one word integer.  The remainder is also a
!  one word integer.

      IMPLICIT NONE

!             Scratch array usage during IMDVIR:   M01 - M03

! .. Intrinsic Functions ..
      INTRINSIC abs, dint, int, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: idiv, irem
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, mda, mdab, mdb, mdr, mkt, modint, mvalp
      INTEGER :: j, jdiv, ka, kl, kltflg, kpt, n1, ndsave, nmval, ntrsav, nv2
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
      EXTERNAL fmwarn, imargs, imdivr, imeq, imi2m, imm2i, imntr, imntri
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMDVIR',1,ma,ma)
      kflag = 0
      ndsave = ndig
      kltflg = 0
      ntrsav = ntrace
      ntrace = 0
      mkt = abs(idiv)
      IF (mkt<mbase) THEN
        m01(0) = ma(0)
        m01(1) = 1
        m01(2) = idiv
        m01(3) = 0
      ELSE IF (mkt<mbase*mbase) THEN
        m01(0) = ma(0)
        m01(1) = 2
        m01(2) = int(mkt/mbase)
        m01(3) = mkt - m01(2)*mbase
        IF (idiv<0) m01(2) = -m01(2)
      ELSE
        CALL imi2m(idiv,m01)
      END IF
      ntrace = ntrsav
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMDVIR'
        CALL imntr(2,ma,ma,1)
        CALL imntri(2,idiv,0)
      END IF
      jdiv = abs(idiv)

!             Check for special cases.

      IF (ma(1)<0) THEN
        mb(1) = munkno
        mb(2) = 1
        mb(3) = 0
        mb(0) = nint(ndg2mx*alogm2)
        irem = iunkno
        kflag = -4
        namest(ncall) = 'IMDVIR'
        CALL fmwarn
        GO TO 70
      END IF
      IF (jdiv==1 .AND. ma(1)/=munkno) THEN
        IF (idiv==1) THEN
          CALL imeq(ma,mb)
          irem = 0
          GO TO 70
        ELSE
          CALL imeq(ma,mb)
          IF (mb(1)/=munkno) mb(2) = -mb(2)
          irem = 0
          GO TO 70
        END IF
      END IF
      IF (ma(1)>ndg2mx .OR. idiv==0) THEN
        kflag = -4
        IF (ma(1)/=munkno) THEN
          namest(ncall) = 'IMDVIR'
          CALL fmwarn
        END IF
        mb(1) = munkno
        mb(2) = 1
        mb(3) = 0
        mb(0) = nint(ndg2mx*alogm2)
        irem = iunkno
        GO TO 70
      END IF
      IF (ma(1)<=2) THEN
        IF (ma(1)<=1) THEN
          mda = ma(2)
        ELSE IF (ma(2)<0) THEN
          mda = ma(2)*mbase - ma(3)
        ELSE
          mda = ma(2)*mbase + ma(3)
        END IF
        mdb = idiv
        mdab = dint(mda/mdb)
        mdr = mda - mdab*mdb
        IF (abs(mdab)<mbase) THEN
          mb(0) = ma(0)
          mb(1) = 1
          IF (mdab==0) mb(1) = 0
          mb(2) = mdab
          mb(3) = 0
        ELSE IF (abs(mdab)<mbase*mbase) THEN
          mb(0) = ma(0)
          mb(1) = 2
          mb(2) = dint(mdab/mbase)
          mb(3) = abs(mdab-mbase*mb(2))
        ELSE
          GO TO 10
        END IF
        irem = int(mdr)
        GO TO 70
      END IF

10    ma2 = ma(2)
      m01(2) = abs(m01(2))
      IF (ma(1)<=m01(1)) THEN
        IF (ma(1)==m01(1) .AND. abs(ma(2))<=m01(2)) THEN
          ma(2) = abs(ma(2))
          IF (imcomp(ma,'EQ',m01)) THEN
            kltflg = 2
          ELSE IF (imcomp(ma,'LT',m01)) THEN
            kltflg = 1
          END IF
          ma(2) = ma2
        END IF
        IF (ma(1)<m01(1) .OR. kltflg>=1) THEN
          IF (kltflg/=2) THEN
            CALL imm2i(ma,irem)
            irem = abs(irem)
            CALL imi2m(0,mb)
          ELSE
            CALL imi2m(1,mb)
            irem = 0
          END IF
          GO TO 60
        END IF
      END IF
      ndig = int(ma(1))
      IF (ndig<2) ndig = 2
      n1 = int(ma(1)) + 1

!             If ABS(IDIV).GE.MXBASE use IMDIVR.

      mvalp = abs(idiv)
      nmval = int(mvalp)
      nv2 = nmval - 1
      IF (abs(idiv)>mxbase .OR. nmval/=abs(idiv) .OR. nv2/=abs(idiv)-1) THEN
        CALL imi2m(idiv,m03)
        CALL imdivr(ma,m03,mb,m03)
        CALL imm2i(m03,irem)
        GO TO 70
      END IF

!             Work with positive numbers.

      ma(2) = abs(ma(2))

!             Find the first significant digit of the quotient.

      mkt = ma(2)
      IF (mkt>=mvalp) THEN
        kpt = 2
        GO TO 30
      END IF
      DO 20 j = 3, n1
        mkt = mkt*mbase + ma(j)
        IF (mkt>=mvalp) THEN
          kpt = j
          GO TO 30
        END IF
20    CONTINUE

      CALL imm2i(ma,irem)
      CALL imi2m(0,mb)
      GO TO 70

!             Do the rest of the division.

30    ka = kpt + 1
      mwa(1) = ma(1) + 2 - kpt
      mwa(2) = int(mkt/mvalp)
      modint = mkt - mwa(2)*mvalp
      IF (ka<=n1) THEN
        kl = 3 - ka

!             (Inner Loop)

        DO 40 j = ka, n1
          mkt = modint*mbase + ma(j)
          mwa(kl+j) = int(mkt/mvalp)
          modint = mkt - mwa(kl+j)*mvalp
40      CONTINUE
      END IF

      mb(0) = ma(0)
      DO 50 j = 1, int(mwa(1)) + 1
        mb(j) = mwa(j)
50    CONTINUE
      irem = int(modint)

60    IF (ma2<0 .AND. idiv>0) THEN
        IF (mb(1)/=munkno) mb(2) = -mb(2)
        irem = -irem
      ELSE IF (ma2>0 .AND. idiv<0) THEN
        IF (mb(1)/=munkno) mb(2) = -mb(2)
      ELSE IF (ma2<0 .AND. idiv<0) THEN
        irem = -irem
      END IF

70    IF (ntrace/=0 .AND. ncall<=lvltrc) THEN
        CALL imntr(1,mb,mb,1)
        CALL imntri(1,irem,0)
      END IF

      ncall = ncall - 1
      ndig = ndsave
      RETURN
    END SUBROUTINE imdvir
    SUBROUTINE imeq(ma,mb)

!  MB = MA

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC int, max
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: j, kdg
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      kdg = max(2,int(ma(1))) + 1
      IF (kdg>lunpck) kdg = 3
      DO 10 j = 0, kdg
        mb(j) = ma(j)
10    CONTINUE
      RETURN
    END SUBROUTINE imeq
    SUBROUTINE imfm2i(ma,mb)

!  MB = INT(MA)

!  Convert from real (FM) format to integer (IM) format.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: ntrsav
! ..
! .. External Subroutines ..
      EXTERNAL fmeq, fmint, fmwarn
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      kflag = 0
      ntrsav = ntrace
      ntrace = 0
      CALL fmeq(ma,mb)
      CALL fmint(mb,mb)
      IF (mb(1)>ndigmx) THEN
        IF (mb(1)<=ndg2mx .OR. ncall<=1) THEN
          mb(0) = nint(ndg2mx*alogm2)
          mb(1) = munkno
          mb(2) = 1
          mb(3) = 0
          kflag = -4
          namest(ncall) = 'IMFM2I'
          CALL fmwarn
        END IF
      END IF
      ntrace = ntrsav
      ncall = ncall - 1

      RETURN
    END SUBROUTINE imfm2i
    SUBROUTINE imform(form,ma,string)

!  Convert an IM number (MA) to a character string base 10 (STRING)
!  using character string FORM format.

!  FORM can be one of these types:  Iw,  Fw.d,  Ew.d,  1PEw.d
!       for positive integers w,d.

      IMPLICIT NONE

!             Scratch array usage during IMFORM:   M01 - M02

! .. Intrinsic Functions ..
      INTRINSIC int
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (*) :: form, string
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmform, imargs
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMFORM',1,ma,ma)
      kflag = 0
      namest(ncall) = 'IMFORM'
      ndsave = ndig
      ndig = int(ma(1))
      IF (ndig<2 .OR. ndig>ndg2mx) ndig = 2
      IF (ma(1)<=1) ma(3) = 0

      CALL fmform(form,ma,string)

      ncall = ncall - 1
      ndig = ndsave
      RETURN
    END SUBROUTINE imform
    SUBROUTINE imfprt(form,ma)

!  Print an IM number (MA) on unit KW using character
!  string FORM format.

!  FORM can be one of these types:  Iw,  Fw.d,  Ew.d,  1PEw.d
!       for positive integers w,d.

      IMPLICIT NONE

!             Scratch array usage during IMFPRT:   M01 - M02

! .. Intrinsic Functions ..
      INTRINSIC int
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (*) :: form
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmfprt, imargs
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMFPRT',1,ma,ma)
      kflag = 0
      namest(ncall) = 'IMFPRT'
      ndsave = ndig
      ndig = int(ma(1))
      IF (ndig<2 .OR. ndig>ndg2mx) ndig = 2
      IF (ma(1)<=1) ma(3) = 0

      CALL fmfprt(form,ma)

      ncall = ncall - 1
      ndig = ndsave
      RETURN
    END SUBROUTINE imfprt
    SUBROUTINE imgcd(ma,mb,mc)

!  MC is returned as the greatest common divisor of MA and MB.

      IMPLICIT NONE

!             Scratch array usage during IMGCD:   M01 - M05

! .. Intrinsic Functions ..
      INTRINSIC abs, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmwarn, imabs, imargs, imdivr, imeq, imi2m, immax, immin, imntr
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMGCD ',2,ma,mb)
      kflag = 0
      ndsave = ndig
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMGCD '
        CALL imntr(2,ma,mb,2)
      END IF

!             Check for special cases.

      IF (ma(1)==munkno .OR. mb(1)==munkno) THEN
        mc(1) = munkno
        mc(2) = 1
        mc(3) = 0
        mc(0) = nint(ndg2mx*alogm2)
        kflag = -4
        GO TO 20
      ELSE IF (mb(2)==0) THEN
        CALL imabs(ma,mc)
        GO TO 20
      ELSE IF (ma(2)==0) THEN
        CALL imabs(mb,mc)
        GO TO 20
      ELSE IF (mb(1)==1 .AND. abs(mb(2))==1) THEN
        CALL imi2m(1,mc)
        GO TO 20
      ELSE IF (ma(1)==1 .AND. abs(ma(2))==1) THEN
        CALL imi2m(1,mc)
        GO TO 20
      ELSE IF (ma(1)>=ndg2mx .OR. mb(1)>=ndg2mx .OR. ma(1)<0 .OR. mb(1)<0) &
          THEN
        mc(1) = munkno
        mc(2) = 1
        mc(3) = 0
        mc(0) = nint(ndg2mx*alogm2)
        kflag = -4
        namest(ncall) = 'IMGCD '
        CALL fmwarn
        GO TO 20
      END IF

      CALL imabs(ma,m05)
      CALL imabs(mb,m04)
      CALL immax(m05,m04,m03)
      CALL immin(m05,m04,m04)
10    CALL imdivr(m03,m04,mc,m05)
      IF (m05(2)/=0) THEN
        CALL imeq(m04,m03)
        CALL imeq(m05,m04)
        GO TO 10
      END IF
      CALL imeq(m04,mc)

      IF (mc(1)==munkno) THEN
        kflag = -4
        namest(ncall) = 'IMGCD '
        CALL fmwarn
      END IF

20    IF (ntrace/=0) CALL imntr(1,mc,mc,1)
      ncall = ncall - 1
      ndig = ndsave
      RETURN
    END SUBROUTINE imgcd
    SUBROUTINE imi2fm(ma,mb)

!  MB = MA

!  Convert from integer (IM) format to real (FM) format.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC int, max, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: kdg
! ..
! .. External Subroutines ..
      EXTERNAL fmequ, imargs
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMI2FM',1,ma,ma)
      kflag = 0
      kdg = max(2,int(ma(1)))
      IF (kdg>ndg2mx) kdg = 2
      IF (ma(1)<=1) ma(3) = 0
      CALL fmequ(ma,mb,kdg,ndig)
      mb(0) = nint(ndg2mx*alogm2)
      ncall = ncall - 1

      RETURN
    END SUBROUTINE imi2fm
    SUBROUTINE imi2m(ival,ma)

!  MA = IVAL

!  Convert a one word integer to IM format.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmim, imntr, imntri
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      kflag = 0
      ndsave = ndig
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMI2M '
        CALL imntri(2,ival,1)

        ndig = 4
        CALL fmim(ival,ma)
        IF (ma(1)>4) THEN
          ndig = ndigmx
          CALL fmim(ival,ma)
        END IF

        CALL imntr(1,ma,ma,1)
      ELSE
        ndig = 4
        CALL fmim(ival,ma)
        IF (ma(1)>4) THEN
          ndig = ndigmx
          CALL fmim(ival,ma)
        END IF
      END IF
      ndig = ndsave
      ncall = ncall - 1
      RETURN
    END SUBROUTINE imi2m
    SUBROUTINE iminp(line,ma,la,lb)

!  Convert an array of characters to multiple precision integer format.

!  LINE is an A1 character array of length LB to be converted
!       to IM format and returned in MA.
!  LA is a pointer telling the routine where in the array to begin
!     the conversion.
!  LB is a pointer to the last character of the field for that number.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC int, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: la, lb
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
      CHARACTER (1) :: line(lb)
! ..
! .. Local Scalars ..
      INTEGER :: kfsave, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fminp, fmint, fmwarn, imntr
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      kflag = 0
      ndsave = ndig
      namest(ncall) = 'IMINP '

      ndig = ndigmx
      CALL fminp(line,ma,la,lb)
      kfsave = kflag
      CALL fmint(ma,ma)
      kflag = kfsave

      IF (ma(1)>ndg2mx .AND. ma(1)<mexpov) THEN
        kflag = -9
        ndig = int(ma(1))
        CALL fmwarn
        ma(0) = nint(ndg2mx*alogm2)
        ma(1) = munkno
        ma(2) = 1
        ma(3) = 0
      END IF

      ndig = ndsave
      IF (ntrace/=0) CALL imntr(1,ma,ma,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE iminp
    SUBROUTINE imm2dp(ma,x)

!  X = MA

!  Convert an IM number to double precision.

!  If KFLAG = -4 is returned for a value of MA that is in the range
!  of the machine's double precision number system, change the
!  definition of DPMAX in routine FMSET to reflect the current machine's
!  range.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dble, int, max
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL (KIND(0.0D0)) :: x
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: kreslt, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmargs, fmmd, fmwarn, imntr, imntrr
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      kflag = 0
      namest(ncall) = 'IMM2DP'
      kreslt = 0
      IF (abs(ma(1))>mexpab) THEN
        CALL fmargs('IMM2DP',1,ma,ma,kreslt)
      END IF
      IF (ntrace/=0) CALL imntr(2,ma,ma,1)
      IF (kreslt/=0) THEN

!             Here no valid result can be returned.  Set X to some
!             value that the user is likely to recognize as wrong.

        x = dble(runkno)
        kflag = -4
        IF (ma(1)/=munkno) CALL fmwarn
        IF (ntrace/=0) CALL imntrr(1,x,1)
        ncall = ncall - 1
        RETURN
      END IF

      ndsave = ndig
      ndig = max(2,int(ma(1)))
      IF (ndig>ndg2mx) ndig = 2
      IF (ma(1)<=1) ma(3) = 0
      CALL fmmd(ma,x)

      IF (ntrace/=0) CALL imntrr(1,x,1)
      ndig = ndsave
      ncall = ncall - 1
      RETURN
    END SUBROUTINE imm2dp
    SUBROUTINE imm2i(ma,ival)

!  IVAL = MA

!  Convert an IM number to a one word integer.

!  KFLAG =  0 is returned if the conversion is exact.
!        = -4 is returned if MA is larger than INTMAX in magnitude.
!             IVAL = IUNKNO is returned as an indication that IVAL
!             could not be computed without integer overflow.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, int
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmm2i, imargs, imntr, imntri
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMM2I ',1,ma,ma)
      ndsave = ndig
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMM2I '
        CALL imntr(2,ma,ma,1)
      END IF

      ndig = int(ma(1))
      IF (ndig<2) ndig = 2
      IF (ndig>ndg2mx) ndig = 2
      IF (ma(1)<=1) ma(3) = 0
      kflag = 0
      CALL fmm2i(ma,ival)

      IF (abs(ntrace)>=1 .AND. ncall<=lvltrc) THEN
        CALL imntri(1,ival,1)
      END IF
      ncall = ncall - 1
      ndig = ndsave
      RETURN
    END SUBROUTINE imm2i
    SUBROUTINE immax(ma,mb,mc)

!  MC = MAX(MA,MB)

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: kwrnsv
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
      EXTERNAL imargs, imeq, imntr
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      kflag = 0
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMMAX ',2,ma,mb)
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMMAX '
        CALL imntr(2,ma,mb,2)
      END IF

      kwrnsv = kwarn
      kwarn = 0
      IF (ma(1)==munkno .OR. mb(1)==munkno) THEN
        mc(1) = munkno
        mc(2) = 1
        mc(0) = nint(ndg2mx*alogm2)
        kflag = -4
      ELSE IF (imcomp(ma,'LT',mb)) THEN
        CALL imeq(mb,mc)
      ELSE
        CALL imeq(ma,mc)
      END IF

      kwarn = kwrnsv
      IF (ntrace/=0) CALL imntr(1,mc,mc,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE immax
    SUBROUTINE immin(ma,mb,mc)

!  MC = MIN(MA,MB)

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: kwrnsv
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
      EXTERNAL imargs, imeq, imntr
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      kflag = 0
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMMIN ',2,ma,mb)
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMMIN '
        CALL imntr(2,ma,mb,2)
      END IF

      kwrnsv = kwarn
      kwarn = 0
      IF (ma(1)==munkno .OR. mb(1)==munkno) THEN
        mc(1) = munkno
        mc(2) = 1
        mc(0) = nint(ndg2mx*alogm2)
        kflag = -4
      ELSE IF (imcomp(ma,'GT',mb)) THEN
        CALL imeq(mb,mc)
      ELSE
        CALL imeq(ma,mc)
      END IF

      kwarn = kwrnsv
      IF (ntrace/=0) CALL imntr(1,mc,mc,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE immin
    SUBROUTINE immod(ma,mb,mc)

!  MC = MOD(MA,MB)

!  Use IMDIVR if both INT(MA/MB) and MOD(MA,MB) are needed.

      IMPLICIT NONE

!             Scratch array usage during IMMOD:   M01 - M03

! .. Intrinsic Functions ..
      INTRINSIC nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmwarn, imargs, imdivr, imntr
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMMOD ',2,ma,mb)
      kflag = 0
      ndsave = ndig
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMMOD '
        CALL imntr(2,ma,mb,2)
      END IF

      IF (ma(1)==munkno .OR. mb(1)==munkno) THEN
        mc(1) = munkno
        mc(2) = 1
        mc(3) = 0
        mc(0) = nint(ndg2mx*alogm2)
        kflag = -4
        GO TO 10
      END IF

      CALL imdivr(ma,mb,m03,mc)

      IF (mc(1)==munkno) THEN
        kflag = -4
        namest(ncall) = 'IMMOD '
        CALL fmwarn
      END IF

10    IF (ntrace/=0) CALL imntr(1,mc,mc,1)
      ncall = ncall - 1
      ndig = ndsave
      RETURN
    END SUBROUTINE immod
    SUBROUTINE immpy(ma,mb,mc)

!  MC = MA * MB

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dint, int, min, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: mdab
      INTEGER :: j, kovfl, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmmpy2, fmwarn, imargs, imeq, imntr
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMMPY ',2,ma,mb)
      kflag = 0
      ndsave = ndig
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMMPY '
        CALL imntr(2,ma,mb,2)
      END IF

      IF (ma(1)<=1) THEN
        IF (mb(1)>1) GO TO 10
        mdab = ma(2)*mb(2)
        IF (abs(mdab)<mbase) THEN
          mc(0) = min(ma(0),mb(0))
          mc(1) = 1
          IF (mdab==0) mc(1) = 0
          mc(2) = mdab
          mc(3) = 0
          GO TO 40
        ELSE IF (abs(mdab)<mbase*mbase) THEN
          mc(0) = min(ma(0),mb(0))
          mc(1) = 2
          mc(2) = dint(mdab/mbase)
          mc(3) = abs(mdab-mbase*mc(2))
          GO TO 40
        END IF
      END IF

!             Check for special cases.

10    kovfl = 0
      IF (ma(1)==mexpov .OR. mb(1)==mexpov) kovfl = 1
      IF (ma(1)<0 .OR. mb(1)<0) THEN
        mc(1) = munkno
        mc(2) = 1
        mc(3) = 0
        mc(0) = nint(ndg2mx*alogm2)
        kflag = -4
        namest(ncall) = 'IMMPY '
        CALL fmwarn
        GO TO 50
      END IF
      IF (mb(1)==1 .AND. mb(2)==1) THEN
        CALL imeq(ma,mc)
        GO TO 40
      ELSE IF (mb(1)==1 .AND. mb(2)==-1) THEN
        CALL imeq(ma,mc)
        IF (mc(1)/=munkno) mc(2) = -mc(2)
        GO TO 40
      ELSE IF (ma(1)==1 .AND. ma(2)==1) THEN
        CALL imeq(mb,mc)
        GO TO 40
      ELSE IF (ma(1)==1 .AND. ma(2)==-1) THEN
        CALL imeq(mb,mc)
        IF (mc(1)/=munkno) mc(2) = -mc(2)
        GO TO 40
      ELSE IF (ma(1)==munkno .OR. mb(1)==munkno) THEN
        kflag = -4
        mc(1) = munkno
        mc(2) = 1
        mc(3) = 0
        mc(0) = nint(ndg2mx*alogm2)
        GO TO 50
      END IF

      ndig = int(ma(1)+mb(1))
      IF (ndig<2) ndig = 2
      IF (ndig>ndg2mx) ndig = ndg2mx
      IF (ma(1)==mexpov .OR. mb(1)==mexpov) ndig = 2

      DO 20 j = int(ma(1)) + 2, ndig + 1
        ma(j) = 0
20    CONTINUE
      DO 30 j = int(mb(1)) + 2, ndig + 1
        mb(j) = 0
30    CONTINUE

      CALL fmmpy2(ma,mb,mc)

      IF (ndig>ndigmx) ndig = 2
40    IF (mc(1)>ndigmx) THEN
        IF (ncall==1 .OR. mc(1)>ndg2mx) THEN
          mc(0) = nint(ndg2mx*alogm2)
          mc(1) = mexpov
          IF (mc(2)>0) THEN
            mc(2) = 1
          ELSE
            mc(2) = -1
          END IF
          mc(3) = 0
          kflag = -5
          namest(ncall) = 'IMMPY '
          IF (kovfl/=1) CALL fmwarn
        END IF
      END IF

50    IF (ntrace/=0) CALL imntr(1,mc,mc,1)
      ncall = ncall - 1
      ndig = ndsave
      RETURN
    END SUBROUTINE immpy
    SUBROUTINE immpy2(ma,mb)

!  Internal multiplication of MA*MB.  The result is returned in MWA.
!  Both MA and MB are positive.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC dint, int, min
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maxmwa, mbj, mbm1, mkt, mmax
      INTEGER :: j, jm1, k, kb, kl, klma, klmb, n1
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      n1 = ndig + 1
      mwa(1) = ma(1) + mb(1)
      mwa(n1+1) = 0

!             The multiplication loop begins here.

!             MAXMWA is an upper bound on the size of values in MWA
!                    divided by (MBASE-1).  It is used to determine
!                    whether to normalize before the next digit is
!                    multiplied.

      mbm1 = mbase - 1
      mmax = intmax - mbase
      mmax = min(dint(maxint/mbm1-mbm1),mmax)
      mbj = mb(2)
      mwa(2) = 0
      klma = int(ma(1))
      DO 10 k = klma + 3, n1
        mwa(k) = 0
10    CONTINUE

!             (Inner Loop)

      DO 20 k = 2, klma + 1
        mwa(k+1) = ma(k)*mbj
20    CONTINUE
      maxmwa = mbj
      klmb = int(mb(1))
      DO 50 j = 3, klmb + 1
        mbj = mb(j)
        IF (mbj/=0) THEN
          maxmwa = maxmwa + mbj
          jm1 = j - 1
          kl = klma + 1

!                       Major (Inner Loop)

          DO 30 k = j + 1, j + klma
            mwa(k) = mwa(k) + ma(k-jm1)*mbj
30        CONTINUE
        END IF

        IF (maxmwa>mmax) THEN
          maxmwa = 0

!                       Here normalization is only required for the
!                       range of digits currently changing in MWA.

          DO 40 kb = jm1 + kl, jm1 + 2, -1
            mkt = int(mwa(kb)/mbase)
            mwa(kb-1) = mwa(kb-1) + mkt
            mwa(kb) = mwa(kb) - mkt*mbase
40        CONTINUE
        END IF
50    CONTINUE

!             Perform the final normalization.  (Inner Loop)

      DO 60 kb = n1, 3, -1
        mkt = int(mwa(kb)/mbase)
        mwa(kb-1) = mwa(kb-1) + mkt
        mwa(kb) = mwa(kb) - mkt*mbase
60    CONTINUE

      RETURN
    END SUBROUTINE immpy2
    SUBROUTINE immpyi(ma,ival,mb)

!  MB = MA * IVAL

!  Multiplication by a one word integer.

      IMPLICIT NONE

!             Scratch array usage during IMMPYI:   M01

! .. Intrinsic Functions ..
      INTRINSIC abs, dble, dint, int, log, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, mcarry, mdab, mkt, mval
      INTEGER :: j, ka, kb, kc, kovfl, kshift, n1, ndsave, nmval, ntrsav, nv2
! ..
! .. External Subroutines ..
      EXTERNAL fmim, fmmpy2, fmwarn, imargs, imeq, imi2m, imntr, imntri
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMMPYI',1,ma,ma)
      kflag = 0
      ndsave = ndig
      ntrsav = ntrace
      ntrace = 0
      ntrace = ntrsav
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMMPYI'
        CALL imntr(2,ma,ma,1)
        CALL imntri(2,ival,0)
      END IF
      ma2 = ma(2)

      IF (ma(1)<=1) THEN
        mdab = ma(2)*ival
        IF (abs(mdab)<mbase) THEN
          mb(0) = ma(0)
          mb(1) = 1
          IF (mdab==0) mb(1) = 0
          mb(2) = mdab
          mb(3) = 0
          GO TO 80
        ELSE IF (abs(mdab)<mbase*mbase) THEN
          mb(0) = ma(0)
          mb(1) = 2
          mb(2) = dint(mdab/mbase)
          mb(3) = abs(mdab-mbase*mb(2))
          GO TO 80
        END IF
      END IF

!             Check for special cases.

      kovfl = 0
      IF (ma(1)==mexpov) kovfl = 1
      IF (ma(1)<0) THEN
        mb(1) = munkno
        mb(2) = 1
        mb(3) = 0
        mb(0) = nint(ndg2mx*alogm2)
        kflag = -4
        namest(ncall) = 'IMMPYI'
        CALL fmwarn
        GO TO 90
      END IF
      IF (ma(1)==munkno) THEN
        kflag = -4
        mb(1) = munkno
        mb(2) = 1
        mb(3) = 0
        mb(0) = nint(ndg2mx*alogm2)
        GO TO 90
      ELSE IF (ival==0) THEN
        CALL imi2m(0,mb)
        GO TO 80
      ELSE IF (ival==1) THEN
        CALL imeq(ma,mb)
        GO TO 80
      ELSE IF (ival==-1) THEN
        CALL imeq(ma,mb)
        IF (mb(1)/=munkno) mb(2) = -mb(2)
        GO TO 80
      ELSE IF (ma(1)==1 .AND. ma(2)==1) THEN
        CALL imi2m(ival,mb)
        GO TO 80
      ELSE IF (ma(1)==1 .AND. ma(2)==-1) THEN
        CALL imi2m(ival,mb)
        IF (mb(1)/=munkno) mb(2) = -mb(2)
        GO TO 80
      ELSE IF (ma(1)==mexpov) THEN
        kflag = -5
        mb(1) = mexpov
        mb(2) = 1
        mb(3) = 0
        mb(0) = nint(ndg2mx*alogm2)
        GO TO 70
      END IF

!             Work with positive numbers.

      ma(2) = abs(ma(2))
      mval = abs(ival)
      nmval = int(mval)
      nv2 = nmval - 1
      ndig = int(ma(1))
      n1 = ndig + 1

!             To leave room for normalization, shift the product
!             to the right KSHIFT places in MWA.

      kshift = int((log(dble(ma(2)+1)*dble(mval)))/dlogmb)

!             If IVAL is too big, use FMMPY.

      IF (kshift>ndig .OR. mval>maxint/mbase .OR. nmval/=abs(ival) .OR. &
          nv2/=abs(ival)-1) THEN
        ma(2) = ma2
        ndig = 4
        CALL fmim(ival,m01)
        ndig = int(ma(1)+m01(1))
        IF (ndig<2) ndig = 2
        IF (ndig>ndg2mx) ndig = ndg2mx
        DO 10 j = int(ma(1)) + 2, ndig + 1
          ma(j) = 0
10      CONTINUE
        IF (ndig>4) CALL fmim(ival,m01)
        CALL fmmpy2(ma,m01,mb)
        GO TO 90
      END IF

      mwa(1) = ma(1) + kshift
      ka = 2 + kshift
      kb = n1 + kshift
      kc = ndig + 5
      DO 20 j = kb, kc
        mwa(j) = 0
20    CONTINUE

      mcarry = 0

!             This is the main multiplication loop.

      DO 30 j = kb, ka, -1
        mkt = ma(j-kshift)*mval + mcarry
        mcarry = int(mkt/mbase)
        mwa(j) = mkt - mcarry*mbase
30    CONTINUE

!             Resolve the final carry.

      DO 40 j = ka - 1, 2, -1
        mkt = int(mcarry/mbase)
        mwa(j) = mcarry - mkt*mbase
        mcarry = mkt
40    CONTINUE

!             Now the first significant digit in the product is in
!             MWA(2) or MWA(3).

      ma(2) = ma2
      mb(0) = ma(0)
      IF (mwa(2)==0) THEN
        mb(1) = mwa(1) - 1
        DO 50 j = 3, kb
          mb(j-1) = mwa(j)
50      CONTINUE
      ELSE
        mb(1) = mwa(1)
        DO 60 j = 2, kb
          mb(j) = mwa(j)
60      CONTINUE
      END IF

!             Put the sign on the result.

70    IF ((ival>0 .AND. ma2<0) .OR. (ival<0 .AND. ma2>0)) mb(2) = -mb(2)

80    IF (mb(1)>ndigmx) THEN
        IF (ncall==1 .OR. mb(1)>ndg2mx) THEN
          mb(0) = nint(ndg2mx*alogm2)
          mb(1) = mexpov
          IF (mb(2)>0) THEN
            mb(2) = 1
          ELSE
            mb(2) = -1
          END IF
          mb(3) = 0
          kflag = -5
          namest(ncall) = 'IMMPYI'
          IF (kovfl/=1) CALL fmwarn
        END IF
      END IF

90    IF (ntrace/=0) CALL imntr(1,mb,mb,1)
      ncall = ncall - 1
      ndig = ndsave
      RETURN
    END SUBROUTINE immpyi
    SUBROUTINE immpym(ma,mb,mc,md)

!  MD = MA * MB mod MC

!  This routine is slightly faster than calling IMMPY and IMMOD
!  separately, and it works for cases where IMMPY would return
!  OVERFLOW.

      IMPLICIT NONE

!             Scratch array usage during IMMPYM:   M01 - M02

! .. Intrinsic Functions ..
      INTRINSIC abs, dble, dint, int, max, min, mod, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, maxmwa, mb2, mbm1, mc1, mc2, mc2p, mcarry, mdab, mdc, &
        mkt, mlmax, mqd
      REAL (KIND(0.0D0)) :: xb, xbase, xbr, xmwa
      INTEGER :: j, jb, jl, k, ka, kb, kl, kltflg, kptmwa, n1, na1, nc1, &
        ndsave, nguard, nl, nmcwds, ntrsav
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
      EXTERNAL fmwarn, imargs, imi2m, immod, immpy2, imntr, imntrj, imprnt, &
        imsub
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMMPYM',2,ma,mb)
      ndsave = ndig
      kflag = 0
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMMPYM'
        CALL imntr(2,ma,mb,2)
        IF (abs(ntrace)>=2 .AND. ncall<=lvltrc) THEN
          IF (ntrace<0) THEN
            ndig = max(2,int(mc(1)))
            IF (ndig>ndg2mx) ndig = 2
            IF (mc(1)<=1) mc(3) = 0
            ntrsav = ntrace
            IF (ntrace<-2) ntrace = -2
            CALL imntrj(mc,ndig)
            ntrace = ntrsav
            ndig = ndsave
          ELSE
            CALL imprnt(mc)
          END IF
        END IF
      END IF

      IF (ma(1)<=1) THEN
        IF (mb(1)>1) GO TO 10
        IF (ma(1)<0 .OR. mb(1)<0) GO TO 10
        mdab = ma(2)*mb(2)
        IF (mc(1)<=2) THEN
          IF (mc(2)==0) GO TO 10
          IF (mc(1)<=1) THEN
            mdc = mc(2)
          ELSE IF (mc(2)<0) THEN
            mdc = mc(2)*mbase - mc(3)
          ELSE
            mdc = mc(2)*mbase + mc(3)
          END IF
          mdab = mod(mdab,mdc)
        END IF
        IF (abs(mdab)<mbase) THEN
          md(0) = min(ma(0),mb(0),mc(0))
          md(1) = 1
          IF (mdab==0) md(1) = 0
          md(2) = mdab
          md(3) = 0
          GO TO 260
        ELSE IF (abs(mdab)<mbase*mbase) THEN
          md(0) = min(ma(0),mb(0),mc(0))
          md(1) = 2
          md(2) = dint(mdab/mbase)
          md(3) = abs(mdab-mbase*md(2))
          GO TO 260
        END IF
      END IF

!             Check for special cases.

10    IF (ma(1)==munkno .OR. mb(1)==munkno .OR. mc(1)==munkno) THEN
        kflag = -4
        md(1) = munkno
        md(2) = 1
        md(3) = 0
        md(0) = nint(ndg2mx*alogm2)
        GO TO 270
      ELSE IF (mc(2)==0 .OR. ma(1)<0 .OR. mb(1)<0 .OR. mc(1)<0) THEN
        kflag = -4
        namest(ncall) = 'IMMPYM'
        CALL fmwarn
        md(1) = munkno
        md(2) = 1
        md(3) = 0
        md(0) = nint(ndg2mx*alogm2)
        GO TO 270
      ELSE IF (ma(2)==0 .OR. mb(2)==0) THEN
        CALL imi2m(0,md)
        GO TO 270
      ELSE IF (mc(1)==1 .AND. abs(mc(2))==1) THEN
        CALL imi2m(0,md)
        GO TO 270
      ELSE IF (mb(1)==1 .AND. mb(2)==1) THEN
        CALL immod(ma,mc,md)
        GO TO 260
      ELSE IF (mb(1)==1 .AND. mb(2)==-1) THEN
        CALL immod(ma,mc,md)
        IF (md(1)/=munkno) md(2) = -md(2)
        GO TO 260
      ELSE IF (ma(1)==1 .AND. ma(2)==1) THEN
        CALL immod(mb,mc,md)
        GO TO 260
      ELSE IF (ma(1)==1 .AND. ma(2)==-1) THEN
        CALL immod(mb,mc,md)
        IF (md(1)/=munkno) md(2) = -md(2)
        GO TO 260
      ELSE IF (ma(1)>ndg2mx .OR. mb(1)>ndg2mx .OR. mc(1)>ndg2mx) THEN
        kflag = -4
        namest(ncall) = 'IMMPYM'
        CALL fmwarn
        md(1) = munkno
        md(2) = 1
        md(3) = 0
        md(0) = nint(ndg2mx*alogm2)
        GO TO 270
      END IF

      ndig = int(ma(1)+mb(1))
      IF (ndig<2) ndig = 2
      IF (ndig>lmwa) ndig = lmwa

!             Save the sign of MA and MB and then work only with
!             positive numbers.

      ma2 = ma(2)
      mb2 = mb(2)
      ma(2) = abs(ma(2))
      mb(2) = abs(mb(2))

      n1 = ndig + 1

!             It is faster if the second argument is the one
!             with fewer digits.

      IF (ma(1)<mb(1)) THEN
        CALL immpy2(mb,ma)
      ELSE
        CALL immpy2(ma,mb)
      END IF

!             Now do the division to find MWA mod MC.

      kltflg = 0
      IF (mwa(2)==0) THEN
        mwa(1) = mwa(1) - 1
      ELSE
        DO 20 j = n1, 2, -1
          mwa(j+1) = mwa(j)
20      CONTINUE
        mwa(2) = 0
      END IF
      mc2 = mc(2)
      kl = int(mc(1))
      IF (kl>lmwa) kl = 2
      DO 30 j = 0, kl + 1
        m01(j) = mc(j)
30    CONTINUE
      m01(2) = abs(m01(2))
      IF (mwa(1)==m01(1) .AND. abs(mwa(3))<=m01(2)) THEN
        DO 40 j = 4, n1
          m02(j-1) = mwa(j)
40      CONTINUE
        m02(2) = abs(mwa(3))
        m02(1) = mwa(1)
        IF (imcomp(m02,'EQ',m01)) THEN
          kltflg = 2
        ELSE IF (imcomp(m02,'LT',m01)) THEN
          kltflg = 1
        END IF
      END IF
      IF (mwa(1)<mc(1) .OR. kltflg>=1) THEN
        IF (kltflg/=2) THEN
          DO 50 j = 3, n1 + 1
            md(j-1) = mwa(j)
50        CONTINUE
          md(1) = mwa(1)
          md(0) = min(ma(0),mb(0),mc(0))
        ELSE
          CALL imi2m(0,md)
        END IF
        GO TO 250
      END IF

      ndig = int(mwa(1))
      IF (ndig<2) ndig = 2

!             NGUARD is the number of guard digits used.

      nguard = 1
      mc2p = abs(mc(2))
      na1 = int(mwa(1)) + 1
      nc1 = int(mc(1)) + 1
      mwa(1) = mwa(1) - mc(1) + 1
      nl = na1 + nguard + 3
      DO 60 j = na1 + 2, nl
        mwa(j) = 0
60    CONTINUE

!             Work only with positive numbers.

      mc1 = mc(1)
      mc2 = mc(2)
      mc(1) = 0
      mc(2) = mc2p

!             NMCWDS is the number of words of MC used to
!             compute the estimated quotient digit MQD.

      nmcwds = 4
      IF (mbase<100) nmcwds = 7

!             XB is an approximation of MC used in
!             estimating the quotient digits.

      xbase = dble(mbase)
      xb = 0
      jl = nmcwds
      IF (jl<=nc1) THEN
        DO 70 j = 2, jl
          xb = xb*xbase + dble(mc(j))
70      CONTINUE
      ELSE
        DO 80 j = 2, jl
          IF (j<=nc1) THEN
            xb = xb*xbase + dble(mc(j))
          ELSE
            xb = xb*xbase
          END IF
80      CONTINUE
      END IF
      IF (jl+1<=nc1) xb = xb + dble(mc(jl+1))/xbase
      xbr = 1.0D0/xb

!             MLMAX determines when to normalize all of MWA.

      mbm1 = mbase - 1
      mlmax = maxint/mbm1
      mkt = intmax - mbase
      mlmax = min(mlmax,mkt)

!             MAXMWA is an upper bound on the size of values in MWA
!             divided by MBASE-1.  It is used to determine whether
!             normalization can be postponed.

      maxmwa = 0

!             KPTMWA points to the next digit in the quotient.

      kptmwa = 2

!             This is the start of the division loop.

!             XMWA is an approximation of the active part of MWA
!             used in estimating quotient digits.

90    kl = kptmwa + nmcwds - 1
      IF (kl<=nl) THEN
        xmwa = ((dble(mwa(kptmwa))*xbase+dble(mwa(kptmwa+1)))*xbase+dble(mwa( &
          kptmwa+2)))*xbase + dble(mwa(kptmwa+3))
        DO 100 j = kptmwa + 4, kl
          xmwa = xmwa*xbase + dble(mwa(j))
100     CONTINUE
      ELSE
        xmwa = dble(mwa(kptmwa))
        DO 110 j = kptmwa + 1, kl
          IF (j<=nl) THEN
            xmwa = xmwa*xbase + dble(mwa(j))
          ELSE
            xmwa = xmwa*xbase
          END IF
110     CONTINUE
      END IF

!             MQD is the estimated quotient digit.

      mqd = dint(xmwa*xbr)
      IF (mqd<0) mqd = mqd - 1

      IF (mqd>0) THEN
        maxmwa = maxmwa + mqd
      ELSE
        maxmwa = maxmwa - mqd
      END IF

!             See if MWA must be normalized.

      ka = kptmwa + 1
      kb = ka + int(mc1) - 1
      IF (maxmwa>=mlmax) THEN
        DO 120 j = kb, ka, -1
          IF (mwa(j)<0) THEN
            mcarry = int((-mwa(j)-1)/mbase) + 1
            mwa(j) = mwa(j) + mcarry*mbase
            mwa(j-1) = mwa(j-1) - mcarry
          ELSE IF (mwa(j)>=mbase) THEN
            mcarry = -int(mwa(j)/mbase)
            mwa(j) = mwa(j) + mcarry*mbase
            mwa(j-1) = mwa(j-1) - mcarry
          END IF
120     CONTINUE
        xmwa = 0
        IF (kl<=nl) THEN
          DO 130 j = kptmwa, kl
            xmwa = xmwa*xbase + dble(mwa(j))
130       CONTINUE
        ELSE
          DO 140 j = kptmwa, kl
            IF (j<=nl) THEN
              xmwa = xmwa*xbase + dble(mwa(j))
            ELSE
              xmwa = xmwa*xbase
            END IF
140       CONTINUE
        END IF
        mqd = dint(xmwa*xbr)
        IF (mqd<0) mqd = mqd - 1
        IF (mqd>0) THEN
          maxmwa = mqd
        ELSE
          maxmwa = -mqd
        END IF
      END IF

!             Subtract MQD*MC from MWA.

      jb = ka - 2
      IF (mqd/=0) THEN

!             Major (Inner Loop)

        DO 150 j = ka, kb
          mwa(j) = mwa(j) - mqd*mc(j-jb)
150     CONTINUE
      END IF

      mwa(ka) = mwa(ka) + mwa(ka-1)*mbase
      mwa(kptmwa) = mqd

      kptmwa = kptmwa + 1
      IF (kptmwa-2<mwa(1)) GO TO 90

!             Final normalization.

      kptmwa = kptmwa - 1
      DO 160 j = kptmwa, 3, -1
        IF (mwa(j)<0) THEN
          mcarry = int((-mwa(j)-1)/mbase) + 1
          mwa(j) = mwa(j) + mcarry*mbase
          mwa(j-1) = mwa(j-1) - mcarry
        ELSE IF (mwa(j)>=mbase) THEN
          mcarry = -int(mwa(j)/mbase)
          mwa(j) = mwa(j) + mcarry*mbase
          mwa(j-1) = mwa(j-1) - mcarry
        END IF
160   CONTINUE

170   DO 180 j = kptmwa + int(mc1), kptmwa + 2, -1
        IF (mwa(j)<0) THEN
          mcarry = int((-mwa(j)-1)/mbase) + 1
          mwa(j) = mwa(j) + mcarry*mbase
          mwa(j-1) = mwa(j-1) - mcarry
        ELSE IF (mwa(j)>=mbase) THEN
          mcarry = -int(mwa(j)/mbase)
          mwa(j) = mwa(j) + mcarry*mbase
          mwa(j-1) = mwa(j-1) - mcarry
        END IF
180   CONTINUE

!             Due to rounding, the remainder may not be between
!             0 and ABS(MC) here.  Correct if necessary.

      IF (mwa(ka)<0) THEN
        DO 190 j = ka, kb
          mwa(j) = mwa(j) + mc(j-jb)
190     CONTINUE
        GO TO 170
      ELSE IF (mwa(ka)>=mbase) THEN
        DO 200 j = ka, kb
          mwa(j) = mwa(j) - mc(j-jb)
200     CONTINUE
        GO TO 170
      END IF

      ma(2) = ma2
      mb(2) = mb2
      mc(1) = mc1
      mc(2) = mc2

      IF (mwa(kptmwa+1)/=0) THEN
        DO 210 j = 1, int(mc1)
          md(j+1) = mwa(kptmwa+j)
210     CONTINUE
        md(1) = mc1
      ELSE
        DO 230 j = 1, int(mc1)
          IF (mwa(kptmwa+j)/=0) THEN
            DO 220 k = j, int(mc1)
              md(k-j+2) = mwa(kptmwa+k)
220         CONTINUE
            md(1) = mc1 + 1 - j
            GO TO 240
          END IF
230     CONTINUE
        md(1) = 0
        md(2) = 0
      END IF
240   IF (md(1)<=1) md(3) = 0
      md(0) = min(ma(0),mb(0),mc(0))

      IF (md(1)>m01(1) .OR. (md(1)==m01(1) .AND. abs(md(2))>=m01(2))) THEN
        IF (imcomp(md,'GE',m01)) CALL imsub(md,m01,md)
      END IF

250   IF (ma2*mb2<0 .AND. mc2>0) THEN
        IF (md(1)/=munkno) md(2) = -md(2)
      ELSE IF (ma2*mb2<0 .AND. mc2<0) THEN
        IF (md(1)/=munkno) md(2) = -md(2)
      END IF

      IF (ndig>ndigmx) ndig = 2
260   IF (md(1)==munkno) THEN
        kflag = -4
        namest(ncall) = 'IMMPYM'
        CALL fmwarn
      END IF

270   IF (ntrace/=0) CALL imntr(1,md,md,1)
      ncall = ncall - 1
      ndig = ndsave
      RETURN
    END SUBROUTINE immpym
    SUBROUTINE imntr(ntr,ma,mb,narg)

!  Print IM numbers in base 10 format.
!  This is used for trace output from the IM routines.

!  NTR =  1 if a result of an IM call is to be printed.
!      =  2 to print input argument(s) to an IM call.

!  MA  -  the IM number to be printed.

!  MB  -  an optional second IM number to be printed.

!  NARG - the number of arguments.  NARG = 1 if only MA is to be
!         printed, and NARG = 2 if both MA and MB are to be printed.

      IMPLICIT NONE

!             Scratch array usage during IMNTR:   M01 - M02

! .. Intrinsic Functions ..
      INTRINSIC abs, int, max
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: narg, ntr
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: ndsave, ntrsav
      CHARACTER (6) :: name
! ..
! .. External Subroutines ..
      EXTERNAL imntrj, imprnt
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (ntrace==0) RETURN
      IF (ncall>lvltrc) RETURN
      IF (ntr==2 .AND. abs(ntrace)==1) RETURN

      IF (ntr==2) THEN
        name = namest(ncall)
        WRITE (kw,90000) name
      ELSE
        name = namest(ncall)
        IF (kflag==0) THEN
          WRITE (kw,90010) name, ncall, int(mbase)
        ELSE
          WRITE (kw,90020) name, ncall, int(mbase), kflag
        END IF
      END IF

      ndsave = ndig
      IF (ntrace<0) THEN
        ndig = max(2,int(ma(1)))
        IF (ndig>ndg2mx) ndig = 2
        IF (ma(1)<=1) ma(3) = 0
        ntrsav = ntrace
        IF (ntrace<-2) ntrace = -2
        CALL imntrj(ma,ndig)
        IF (narg==2) THEN
          ndig = max(2,int(mb(1)))
          IF (ndig>ndg2mx) ndig = 2
          IF (mb(1)<=1) mb(3) = 0
          CALL imntrj(mb,ndig)
        END IF
        ntrace = ntrsav
      END IF

      IF (ntrace>0) THEN
        CALL imprnt(ma)
        IF (narg==2) CALL imprnt(mb)
      END IF

      ndig = ndsave
      RETURN
90000 FORMAT (' Input to ',A6)
90010 FORMAT (' ',A6,15X,'Call level =',I2,5X,'MBASE =',I10)
90020 FORMAT (' ',A6,6X,'Call level =',I2,4X,'MBASE =',I10,4X,'KFLAG =',I3)
    END SUBROUTINE imntr
    SUBROUTINE imntri(ntr,n,knam)

!  Internal routine for trace output of integer variables.

!  NTR = 1 for output values
!        2 for input values

!  N     Integer to be printed.

!  KNAM  is positive if the routine name is to be printed.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, int
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: knam, n, ntr
! ..
! .. Local Scalars ..
      CHARACTER (6) :: name
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (ntrace==0) RETURN
      IF (ncall>lvltrc) RETURN
      IF (ntr==2 .AND. abs(ntrace)==1) RETURN

      IF (ntr==2 .AND. knam>0) THEN
        name = namest(ncall)
        WRITE (kw,90000) name
      END IF
      IF (ntr==1 .AND. knam>0) THEN
        name = namest(ncall)
        IF (kflag==0) THEN
          WRITE (kw,90010) name, ncall, int(mbase)
        ELSE
          WRITE (kw,90020) name, ncall, int(mbase), kflag
        END IF
      END IF

      WRITE (kw,90030) n

      RETURN
90000 FORMAT (' Input to ',A6)
90010 FORMAT (' ',A6,15X,'Call level =',I2,5X,'MBASE =',I10)
90020 FORMAT (' ',A6,6X,'Call level =',I2,4X,'MBASE =',I10,4X,'KFLAG =',I3)
90030 FORMAT (1X,I18)
    END SUBROUTINE imntri
    SUBROUTINE imntrj(ma,nd)

!  Print trace output in internal base MBASE format.  The number to
!  be printed is in MA.

!  ND is the number of base MBASE digits to be printed.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC dble, int, log10
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: nd
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: j, l, n, n1
      CHARACTER (50) :: form
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      n1 = nd + 1

      l = int(log10(dble(mbase-1))) + 2
      n = (kswide-23)/l
      IF (n>10) n = 5*(n/5)
      IF (nd<=n) THEN
        WRITE (form,90000) l + 2, n - 1, l
      ELSE
        WRITE (form,90010) l + 2, n - 1, l, n, l
      END IF
      WRITE (kw,form) (int(ma(j)),j=1,n1)

      RETURN
90000 FORMAT (' (1X,I19,I',I2,',',I3,'I',I2,') ')
90010 FORMAT (' (1X,I19,I',I2,',',I3,'I',I2,'/(22X,',I3,'I',I2,')) ')
    END SUBROUTINE imntrj
    SUBROUTINE imntrr(ntr,x,knam)

!  Internal routine for trace output of real variables.

!  NTR - 1 for output values
!        2 for input values

!  X   - Double precision value to be printed if NX.EQ.1

!  KNAM - Positive if the routine name is to be printed.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, int
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL (KIND(0.0D0)) :: x
      INTEGER :: knam, ntr
! ..
! .. Local Scalars ..
      CHARACTER (6) :: name
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (ntrace==0) RETURN
      IF (ncall>lvltrc) RETURN
      IF (ntr==2 .AND. abs(ntrace)==1) RETURN

      IF (ntr==2 .AND. knam>0) THEN
        name = namest(ncall)
        WRITE (kw,90000) name
      END IF
      IF (ntr==1 .AND. knam>0) THEN
        name = namest(ncall)
        IF (kflag==0) THEN
          WRITE (kw,90010) name, ncall, int(mbase)
        ELSE
          WRITE (kw,90020) name, ncall, int(mbase), kflag
        END IF
      END IF

      WRITE (kw,90030) x

      RETURN
90000 FORMAT (' Input to ',A6)
90010 FORMAT (' ',A6,15X,'Call level =',I2,5X,'MBASE =',I10)
90020 FORMAT (' ',A6,6X,'Call level =',I2,4X,'MBASE =',I10,4X,'KFLAG =',I3)
90030 FORMAT (1X,D30.20)
    END SUBROUTINE imntrr
    SUBROUTINE imout(ma,line,lb)

!  Convert an integer multiple precision number to a character array
!  for output.

!  MA   is an IM number to be converted to an A1 character
!       array in base 10 format
!  LINE is the CHARACTER*1 array in which the result is returned.
!  LB   is the length of LINE.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC int, max
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: lb
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
      CHARACTER (1) :: line(lb)
! ..
! .. Local Scalars ..
      INTEGER :: jf1sav, jf2sav, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmout, imargs
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMOUT ',1,ma,ma)
      kflag = 0
      ndsave = ndig
      namest(ncall) = 'IMOUT '

      ndsave = ndig
      jf1sav = jform1
      jf2sav = jform2
      jform1 = 2
      jform2 = 0
      ndig = max(2,int(ma(1)))
      IF (ndig>ndg2mx) ndig = 2
      IF (ma(1)<=1) ma(3) = 0
      CALL fmout(ma,line,lb)

      ndig = ndsave
      jform1 = jf1sav
      jform2 = jf2sav
      ncall = ncall - 1
      RETURN
    END SUBROUTINE imout
    SUBROUTINE impack(ma,mp)

!  MA is packed two base NDIG digits per word and returned in MP.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, int, mod
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mp(0:lpack)
! ..
! .. Local Scalars ..
      INTEGER :: j, kma1, kp
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      kma1 = int(ma(1))
      IF (kma1<=2 .OR. kma1>ndg2mx) kma1 = 2
      kp = 2
      mp(0) = ma(0)
      mp(1) = ma(1)
      mp(2) = abs(ma(2))*mbase + ma(3)
      IF (ma(2)<0) mp(2) = -mp(2)
      IF (kma1>=4) THEN
        DO 10 j = 4, kma1, 2
          kp = kp + 1
          mp(kp) = ma(j)*mbase + ma(j+1)
10      CONTINUE
      END IF
      IF (mod(kma1,2)==1) mp(kp+1) = ma(kma1+1)*mbase
      RETURN
    END SUBROUTINE impack
    SUBROUTINE impmod(ma,mb,mc,md)

!  MD = MOD(MA**MB,MC)

!  The binary multiplication method used requires an average of
!  1.5 * LOG2(MB) operations.

      IMPLICIT NONE

!             Scratch array usage during IMPMOD:   M01 - M06

! .. Intrinsic Functions ..
      INTRINSIC abs, int, max, min, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck), md(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: macca, maccb, mb2
      INTEGER :: irem, kwrnsv, ndsave, ntrsav
! ..
! .. External Subroutines ..
      EXTERNAL fmwarn, imabs, imargs, imdivr, imdvir, imeq, imi2m, immod, &
        immpym, imntr, imntrj, imprnt
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      kflag = 0
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMPMOD',2,ma,mb)
      IF (kdebug==1) CALL imargs('IMPMOD',1,mc,mc)
      ndsave = ndig
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMPMOD'
        CALL imntr(2,ma,mb,2)
        IF (abs(ntrace)>=2 .AND. ncall<=lvltrc) THEN
          IF (ntrace<0) THEN
            ndig = max(2,int(mc(1)))
            IF (ndig>ndg2mx) ndig = 2
            IF (mc(1)<=1) mc(3) = 0
            ntrsav = ntrace
            IF (ntrace<-2) ntrace = -2
            CALL imntrj(mc,ndig)
            ntrace = ntrsav
            ndig = ndsave
          ELSE
            CALL imprnt(mc)
          END IF
        END IF
      END IF
      mb2 = mb(2)
      macca = ma(0)
      maccb = mb(0)

!             Check for special cases.

      IF (ma(1)==munkno .OR. mb(1)==munkno .OR. mc(1)==munkno .OR. &
          ma(1)==mexpov .OR. mb(1)==mexpov .OR. mc(1)==mexpov .OR. &
          ma(1)<0 .OR. mb(1)<0 .OR. mc(1)<0 .OR. (mb(2)<=0 .AND. ma( &
          2)==0) .OR. mc(2)==0) THEN
        kflag = -4
        IF (ma(1)/=munkno .AND. mb(1)/=munkno .AND. mc(1)/=munkno) THEN
          namest(ncall) = 'IMPMOD'
          CALL fmwarn
        END IF
        md(0) = nint(ndg2mx*alogm2)
        md(1) = munkno
        md(2) = 1
        md(3) = 0
        IF (ntrace/=0) CALL imntr(1,md,md,1)
        ncall = ncall - 1
        RETURN
      END IF

      IF (mb2==0) THEN
        CALL imi2m(1,md)
        IF (ntrace/=0) CALL imntr(1,md,md,1)
        ncall = ncall - 1
        RETURN
      END IF

      IF (mb(1)==1 .AND. abs(mb2)==1) THEN
        kwrnsv = kwarn
        kwarn = 0
        IF (mb2==1) THEN
          CALL immod(ma,mc,md)
        ELSE
          CALL imi2m(1,m05)
          CALL imdivr(m05,ma,m04,m06)
          CALL immod(m04,mc,md)
        END IF
        IF (ntrace/=0) CALL imntr(1,md,md,1)
        ncall = ncall - 1
        kwarn = kwrnsv
        RETURN
      END IF

      IF (ma(2)==0) THEN
        CALL imi2m(0,md)
        IF (ntrace/=0) CALL imntr(1,md,md,1)
        ncall = ncall - 1
        RETURN
      END IF

!             Initialize.

      kwrnsv = kwarn
      kwarn = 0
      CALL imabs(mb,m06)
      CALL imdivr(ma,mc,m04,m05)
      CALL imeq(mc,m04)
      CALL imdvir(m06,2,md,irem)
      IF (irem==0) THEN
        CALL imi2m(1,md)
      ELSE
        CALL imeq(m05,md)
      END IF
      CALL imdvir(m06,2,m06,irem)

!             This is the multiplication loop.

10    CALL imdvir(m06,2,m06,irem)
      CALL immpym(m05,m05,m04,m05)
      IF (irem==1) CALL immpym(m05,md,m04,md)
      IF (m06(2)>0 .AND. md(2)/=0) GO TO 10

      IF (mb2<0) THEN
        CALL imi2m(1,m05)
        CALL imdivr(m05,md,md,m06)
      END IF
      kwarn = kwrnsv
      md(0) = min(macca,maccb)
      IF (kflag<0) THEN
        namest(ncall) = 'IMPMOD'
        CALL fmwarn
      END IF
      IF (ntrace/=0) CALL imntr(1,md,md,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE impmod
    SUBROUTINE imprnt(ma)

!  Print MA in base 10 format.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC int, max
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: jf1sav, jf2sav, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmprnt
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ndsave = ndig
      jf1sav = jform1
      jf2sav = jform2
      jform1 = 2
      jform2 = 0
      ndig = max(2,int(ma(1)))
      IF (ma(1)<=1) ma(3) = 0
      IF (ndig>ndg2mx) ndig = 2
      CALL fmprnt(ma)
      jform1 = jf1sav
      jform2 = jf2sav
      ndig = ndsave
      RETURN
    END SUBROUTINE imprnt
    SUBROUTINE impwr(ma,mb,mc)

!  MC = MA ** MB

!  The binary multiplication method used requires an average of
!  1.5 * LOG2(MB) multiplications.

      IMPLICIT NONE

!             Scratch array usage during IMPWR:   M01 - M06

! .. Intrinsic Functions ..
      INTRINSIC abs, min, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, macca, maccb, mb2
      INTEGER :: irem, iremb, jsign, kwrnsv
! ..
! .. External Subroutines ..
      EXTERNAL fmwarn, imabs, imargs, imdivr, imdvir, imeq, imi2m, immpy, &
        imntr, imsqr
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mesav(0:lunpck), mlbsav(0:lunpck), &
        mln1(0:lunpck), mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), &
        mpisav(0:lunpck), mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      kflag = 0
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMPWR ',2,ma,mb)
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMPWR '
        CALL imntr(2,ma,mb,2)
      END IF
      ma2 = ma(2)
      mb2 = mb(2)
      macca = ma(0)
      maccb = mb(0)
      kwrnsv = kwarn

!             Check for special cases.

      IF (ma(1)==munkno .OR. mb(1)==munkno .OR. ma(1)<0 .OR. mb(1)<0 .OR. (mb( &
          2)<=0 .AND. ma(2)==0)) THEN
        kflag = -4
        IF (ma(1)/=munkno .AND. mb(1)/=munkno) THEN
          kwarn = kwrnsv
          namest(ncall) = 'IMPWR '
          CALL fmwarn
        END IF
        mc(0) = nint(ndg2mx*alogm2)
        mc(1) = munkno
        mc(2) = 1
        mc(3) = 0
        GO TO 30
      END IF

      IF (mb2==0) THEN
        CALL imi2m(1,mc)
        GO TO 30
      END IF

      IF (ma(1)==1 .AND. abs(ma2)==1) THEN
        kwarn = 0
        IF (ma2==1) THEN
          CALL imi2m(1,mc)
        ELSE
          CALL imi2m(2,m05)
          CALL imdivr(mb,m05,m05,m06)
          IF (m06(1)==munkno) THEN
            mc(0) = nint(ndg2mx*alogm2)
            mc(1) = munkno
            mc(2) = 1
            mc(3) = 0
            kflag = -4
            kwarn = kwrnsv
            namest(ncall) = 'IMPWR '
            CALL fmwarn
          ELSE IF (m06(2)==0) THEN
            CALL imi2m(1,mc)
          ELSE
            CALL imi2m(-1,mc)
          END IF
        END IF
        GO TO 30
      END IF

      IF (mb(1)==1 .AND. abs(mb2)==1) THEN
        kwarn = 0
        IF (mb2==1) THEN
          CALL imeq(ma,mc)
        ELSE
          CALL imi2m(1,m05)
          CALL imdivr(m05,ma,mc,m06)
        END IF
        GO TO 30
      END IF

      IF (ma(2)==0) THEN
        CALL imi2m(0,mc)
        GO TO 30
      END IF

      IF (mb(1)==mexpov) THEN
        IF (mb2<0) THEN
          CALL imi2m(0,mc)
        ELSE IF (ma2>0) THEN
          mc(0) = nint(ndg2mx*alogm2)
          mc(1) = mexpov
          mc(2) = 1
          mc(3) = 0
          kflag = -5
        ELSE
          mc(0) = nint(ndg2mx*alogm2)
          mc(1) = munkno
          mc(2) = 1
          mc(3) = 0
          kflag = -4
          kwarn = kwrnsv
          namest(ncall) = 'IMPWR '
          CALL fmwarn
        END IF
        GO TO 30
      END IF

      IF (ma(1)==mexpov) THEN
        jsign = 1
        IF (ma(2)<0) jsign = -1
        IF (mb2>0) THEN
          CALL imdvir(mb,2,mc,irem)
          mc(0) = nint(ndg2mx*alogm2)
          mc(1) = mexpov
          mc(2) = jsign**irem
          mc(3) = 0
          kflag = -5
        ELSE
          CALL imi2m(0,mc)
        END IF
        GO TO 30
      END IF

!             Initialize.

      kwarn = 0
      CALL imabs(mb,m06)

      CALL imeq(ma,m05)

      CALL imdvir(mb,2,mc,iremb)
      IF (iremb==0) THEN
        CALL imi2m(1,mc)
      ELSE
        CALL imeq(m05,mc)
      END IF
      CALL imdvir(m06,2,m06,irem)

!             This is the multiplication loop.

10    CALL imdvir(m06,2,m06,irem)
      CALL imsqr(m05,m05)
      IF (irem==1) CALL immpy(m05,mc,mc)
      IF (m05(1)==mexpov) THEN
        CALL imeq(m05,mc)
        IF (ma2<0 .AND. iremb==1) mc(2) = -1
        GO TO 20
      END IF
      IF (m06(2)>0) GO TO 10

20    IF (mb2<0) THEN
        CALL imi2m(1,m05)
        CALL imdivr(m05,mc,mc,m06)
      END IF

      mc(0) = min(macca,maccb)
      IF (mc(1)>ndigmx) THEN
        IF (ncall==1 .OR. mc(1)>ndg2mx) THEN
          mc(0) = nint(ndg2mx*alogm2)
          mc(1) = mexpov
          IF (mc(2)>0) THEN
            mc(2) = 1
          ELSE
            mc(2) = -1
          END IF
          mc(3) = 0
          kflag = -5
          kwarn = kwrnsv
          namest(ncall) = 'IMPWR '
          CALL fmwarn
        END IF
      END IF

30    kwarn = kwrnsv
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMPWR '
        CALL imntr(1,mc,mc,1)
      END IF
      ncall = ncall - 1
      RETURN
    END SUBROUTINE impwr
    SUBROUTINE imread(kread,ma)

!  Read MA on unit KREAD.  Multi-line numbers will have '&' as the
!  last nonblank character on all but the last line.  Only one
!  number is allowed on the line(s).

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kread
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: kwrnsv, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmnint, fmread
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      ndsave = ndig
      ndig = ndigmx
      CALL fmread(kread,ma)
      kwrnsv = kwarn
      kwarn = 0
      CALL fmnint(ma,ma)
      kwarn = kwrnsv
      ndig = ndsave
      ncall = ncall - 1
      RETURN
    END SUBROUTINE imread
    SUBROUTINE imsign(ma,mb,mc)

!  MC = SIGN(MA,MB)

!  MC is set to ABS(MA) if MB is positive or zero,
!     or -ABS(MA) if MB is negative.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, int, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: kwrnsv, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmwarn, imargs, imeq, imntr
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      kflag = 0
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMSIGN',2,ma,mb)
      ndsave = ndig
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMSIGN'
        CALL imntr(2,ma,mb,2)
      END IF

      ndig = int(ma(1))
      IF (ndig<2) ndig = 2
      IF (ndig>ndg2mx) ndig = 2
      IF (ma(1)<=1) ma(3) = 0
      kwrnsv = kwarn
      kwarn = 0
      IF (ma(1)==munkno .OR. mb(1)==munkno) THEN
        mc(1) = munkno
        mc(2) = 1
        mc(3) = 0
        mc(0) = nint(ndg2mx*alogm2)
        kflag = -4
      ELSE IF (ma(1)<0 .OR. mb(1)<0) THEN
        mc(1) = munkno
        mc(2) = 1
        mc(3) = 0
        mc(0) = nint(ndg2mx*alogm2)
        kflag = -4
        namest(ncall) = 'IMSIGN'
        CALL fmwarn
      ELSE IF (mb(2)>=0) THEN
        CALL imeq(ma,mc)
        mc(2) = abs(mc(2))
      ELSE
        CALL imeq(ma,mc)
        mc(2) = -abs(mc(2))
      END IF

      kwarn = kwrnsv
      IF (ntrace/=0) CALL imntr(1,mc,mc,1)
      ncall = ncall - 1
      ndig = ndsave
      RETURN
    END SUBROUTINE imsign
    SUBROUTINE imsqr(ma,mb)

!  MB = MA * MA

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dint, int, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: mdab
      INTEGER :: j, kovfl, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmsqr2, fmwarn, imargs, imi2m, imntr
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMSQR ',1,ma,ma)
      kflag = 0
      ndsave = ndig
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMSQR '
        CALL imntr(2,ma,ma,1)
      END IF

      IF (ma(1)<=1) THEN
        IF (ma(1)<0) GO TO 10
        mdab = ma(2)*ma(2)
        IF (abs(mdab)<mbase) THEN
          mb(0) = ma(0)
          mb(1) = 1
          IF (mdab==0) mb(1) = 0
          mb(2) = mdab
          mb(3) = 0
          GO TO 30
        ELSE IF (abs(mdab)<mbase*mbase) THEN
          mb(0) = ma(0)
          mb(1) = 2
          mb(2) = dint(mdab/mbase)
          mb(3) = abs(mdab-mbase*mb(2))
          GO TO 30
        END IF
      END IF

!             Check for special cases.

10    kovfl = 0
      IF (ma(1)==mexpov) kovfl = 1
      IF (ma(1)==1 .AND. abs(ma(2))==1) THEN
        CALL imi2m(1,mb)
        GO TO 30
      ELSE IF (ma(1)==munkno) THEN
        kflag = -4
        mb(1) = munkno
        mb(2) = 1
        mb(3) = 0
        mb(0) = nint(ndg2mx*alogm2)
        GO TO 40
      ELSE IF (ma(1)<0) THEN
        kflag = -4
        mb(1) = munkno
        mb(2) = 1
        mb(3) = 0
        mb(0) = nint(ndg2mx*alogm2)
        namest(ncall) = 'IMSQR '
        CALL fmwarn
        GO TO 40
      END IF

      ndig = int(ma(1)+ma(1))
      IF (ndig<2) ndig = 2
      IF (ndig>ndg2mx) ndig = ndg2mx
      IF (ma(1)==mexpov .OR. mb(1)==mexpov) ndig = 2

      DO 20 j = int(ma(1)) + 2, ndig + 1
        ma(j) = 0
20    CONTINUE

      CALL fmsqr2(ma,mb)

      IF (ndig>ndigmx) ndig = 2
30    IF (mb(1)>ndigmx) THEN
        IF (ncall==1 .OR. mb(1)>ndg2mx) THEN
          mb(0) = nint(ndg2mx*alogm2)
          mb(1) = mexpov
          IF (mb(2)>0) THEN
            mb(2) = 1
          ELSE
            mb(2) = -1
          END IF
          mb(3) = 0
          kflag = -5
          namest(ncall) = 'IMSQR '
          IF (kovfl/=1) CALL fmwarn
        END IF
      END IF

40    IF (ntrace/=0) CALL imntr(1,mb,mb,1)
      ncall = ncall - 1
      ndig = ndsave
      RETURN
    END SUBROUTINE imsqr
    SUBROUTINE imst2m(string,ma)

!  MA = STRING

!  Convert a character string to IM format.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC len
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (*) :: string
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: j, lb
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, iminp
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      IF (mblogs/=mbase) CALL fmcons
      ncall = ncall + 1
      namest(ncall) = 'IMST2M'
      lb = len(string)

      DO 10 j = 1, lb
        cmbuff(j) = string(j:j)
10    CONTINUE

      CALL iminp(cmbuff,ma,1,lb)

      ncall = ncall - 1
      RETURN
    END SUBROUTINE imst2m
    SUBROUTINE imsub(ma,mb,mc)

!  MC = MA - MB

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dint, int, min, nint, sign
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mb(0:lunpck), mc(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: mda, mdab, mdb
      INTEGER :: j, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmadd2, fmwarn, imargs, imntr
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMSUB ',2,ma,mb)
      kflag = 0
      ndsave = ndig
      IF (ntrace/=0) THEN
        namest(ncall) = 'IMSUB '
        CALL imntr(2,ma,mb,2)
      END IF

      IF (ma(1)<=2) THEN
        IF (mb(1)>2 .OR. ma(1)<0 .OR. mb(1)<0) GO TO 10
        IF (ma(1)<=1) THEN
          mda = ma(2)
        ELSE IF (ma(2)<0) THEN
          mda = ma(2)*mbase - ma(3)
        ELSE
          mda = ma(2)*mbase + ma(3)
        END IF
        IF (mb(1)<=1) THEN
          mdb = mb(2)
        ELSE IF (mb(2)<0) THEN
          mdb = mb(2)*mbase - mb(3)
        ELSE
          mdb = mb(2)*mbase + mb(3)
        END IF
        mdab = mda - mdb
        IF (abs(mdab)<mbase) THEN
          mc(0) = min(ma(0),mb(0))
          mc(1) = 1
          IF (mdab==0) mc(1) = 0
          mc(2) = mdab
          mc(3) = 0
          IF (mda==0 .OR. mdb==0) kflag = 1
          GO TO 40
        ELSE IF (abs(mdab)<mbase*mbase) THEN
          mc(0) = min(ma(0),mb(0))
          mc(1) = 2
          mc(2) = dint(mdab/mbase)
          mc(3) = abs(mdab-mbase*mc(2))
          IF (mda==0 .OR. mdb==0) kflag = 1
          GO TO 40
        END IF
      END IF

!             Check for special cases.

10    IF (ma(1)>ndg2mx .OR. mb(1)>ndg2mx .OR. ma(1)<0 .OR. mb(1)<0) THEN
        IF (ma(1)==munkno .OR. mb(1)==munkno) THEN
          mc(1) = munkno
          mc(2) = 1
          mc(3) = 0
          mc(0) = nint(ndg2mx*alogm2)
          kflag = -4
          GO TO 50
        END IF
        IF (ma(1)==mexpov) THEN
          mda = 1
          IF ((sign(mda,ma(2))==sign(mda,-mb(2))) .OR. (mb(2)==0)) THEN
            mc(0) = ma(0)
            mc(1) = ma(1)
            mc(2) = ma(2)
            mc(3) = ma(3)
            kflag = -5
            GO TO 50
          ELSE
            mc(0) = nint(ndg2mx*alogm2)
            mc(1) = munkno
            mc(2) = 1
            mc(3) = 0
            kflag = -4
            namest(ncall) = 'IMSUB '
            CALL fmwarn
            GO TO 50
          END IF
        END IF
        IF (mb(1)==mexpov) THEN
          mda = 1
          IF ((sign(mda,-mb(2))==sign(mda,ma(2))) .OR. (ma(2)==0)) THEN
            mc(0) = mb(0)
            mc(1) = mb(1)
            mc(2) = -mb(2)
            mc(3) = mb(3)
            kflag = -5
            GO TO 50
          ELSE
            mc(0) = nint(ndg2mx*alogm2)
            mc(1) = munkno
            mc(2) = 1
            mc(3) = 0
            kflag = -4
            namest(ncall) = 'IMSUB '
            CALL fmwarn
            GO TO 50
          END IF
        END IF
        mc(1) = munkno
        mc(2) = 1
        mc(3) = 0
        mc(0) = nint(ndg2mx*alogm2)
        kflag = -4
        namest(ncall) = 'IMSUB '
        CALL fmwarn
        GO TO 50
      END IF

      IF (ma(1)>mb(1)) THEN
        ndig = int(ma(1)) + 1
        IF (ndig<2 .OR. ndig>ndg2mx) ndig = 2
        ma(ndig+1) = 0
        DO 20 j = int(mb(1)) + 2, ndig + 1
          mb(j) = 0
20      CONTINUE
      ELSE
        ndig = int(mb(1)) + 1
        IF (ndig<2 .OR. ndig>ndg2mx) ndig = 2
        mb(ndig+1) = 0
        DO 30 j = int(ma(1)) + 2, ndig + 1
          ma(j) = 0
30      CONTINUE
      END IF

!             FMADD2 will negate MB and add.

      ksub = 1
      CALL fmadd2(ma,mb,mc)
      ksub = 0

40    IF (mc(1)>ndigmx) THEN
        IF (ncall==1 .OR. mc(1)>ndg2mx) THEN
          mc(0) = nint(ndg2mx*alogm2)
          mc(1) = mexpov
          IF (mc(2)>0) THEN
            mc(2) = 1
          ELSE
            mc(2) = -1
          END IF
          mc(3) = 0
          kflag = -5
          namest(ncall) = 'IMSUB '
          CALL fmwarn
        END IF
      END IF

50    IF (ntrace/=0) CALL imntr(1,mc,mc,1)
      ncall = ncall - 1
      ndig = ndsave
      RETURN
    END SUBROUTINE imsub
    SUBROUTINE imunpk(mp,ma)

!  MP is unpacked and the value returned in MA.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dint, int, mod
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck), mp(0:lpack)
! ..
! .. Local Scalars ..
      INTEGER :: j, kma1, kp
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ntrace
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      kma1 = int(mp(1))
      IF (kma1<=2 .OR. kma1>ndg2mx) kma1 = 2
      kp = 2
      ma(0) = mp(0)
      ma(1) = mp(1)
      ma(2) = dint(abs(mp(2))/mbase)
      ma(3) = abs(mp(2)) - ma(2)*mbase
      IF (mp(2)<0) ma(2) = -ma(2)
      IF (kma1>=4) THEN
        DO 10 j = 4, kma1, 2
          kp = kp + 1
          ma(j) = dint(mp(kp)/mbase)
          ma(j+1) = mp(kp) - ma(j)*mbase
10      CONTINUE
      END IF
      IF (mod(kma1,2)==1) ma(kma1+1) = dint(mp(kp+1)/mbase)
      RETURN
    END SUBROUTINE imunpk
    SUBROUTINE imwrit(kwrite,ma)

!  Write MA on unit KWRITE.  Multi-line numbers will have '&' as the
!  last nonblank character on all but the last line.  These numbers can
!  then be read easily using IMREAD.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC int, log10, max, min, mod, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kwrite
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: j, k, ksave, l, last, lb, nd, ndsave, nexp
! ..
! .. External Subroutines ..
      EXTERNAL imargs, imout
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, kaccsw, kdebug, keswch, &
        kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, ncall, ndg2mx, &
        ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmsave/ndigpi, ndige, ndiglb, ndigli, mbspi, mbse, mbslb, mbsli, &
        mpisav, mesav, mlbsav, mln1, mln2, mln3, mln4, mblogs, mexpab, alogmb, &
        alogm2, alogmx, alogmt, dlogmb, dlogtn, dlogtw, dlogtp, dlogpi, dppi, &
        dpeps, dlogeb, khasht, khashv, ngrd21, ngrd52, ngrd22
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
      ncall = ncall + 1
      IF (kdebug==1) CALL imargs('IMWRIT',1,ma,ma)
      namest(ncall) = 'IMWRIT'
      ndsave = ndig
      ndig = max(2,int(ma(1)))
      IF (ndig>ndg2mx) ndig = 2

      ksave = kflag
      nd = int(real(ndig)*log10(real(mbase))) + 1
      IF (nd<2) nd = 2
      nexp = int(2.0*log10(real(mxbase))) + 6
      lb = min(nd+nexp,lmbuff)

      CALL imout(ma,cmbuff,lb)

      kflag = ksave
      ndig = ndsave
      last = lb + 1
      DO 10 j = 1, lb
        IF (cmbuff(last-j)/=' ' .OR. j==lb) THEN
          l = last - j
          IF (mod(l,73)/=0) THEN
            WRITE (kwrite,90000) (cmbuff(k),k=1,l)
          ELSE
            IF (l>73) WRITE (kwrite,90000) (cmbuff(k),k=1,l-73)
            WRITE (kwrite,90010) (cmbuff(k),k=l-72,l)
          END IF
          ncall = ncall - 1
          RETURN
        END IF
10    CONTINUE
      ncall = ncall - 1
      RETURN
90000 FORMAT (4X,73A1,' &')
90010 FORMAT (4X,73A1)
    END SUBROUTINE imwrit

!  These versions of the IM routines use packed IM numbers.

    SUBROUTINE ipabs(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL imabs, impack, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imabs(mx,my)
      CALL impack(my,mb)
      RETURN
    END SUBROUTINE ipabs
    SUBROUTINE ipadd(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL imadd, impack, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imunpk(mb,my)
      CALL imadd(mx,my,mx)
      CALL impack(mx,mc)
      RETURN
    END SUBROUTINE ipadd
    SUBROUTINE ipbig(ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL imbig, impack
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imbig(my)
      CALL impack(my,ma)
      RETURN
    END SUBROUTINE ipbig
    FUNCTION ipcomp(ma,lrel,mb)
      IMPLICIT NONE
! .. Function Return Value ..
      LOGICAL :: ipcomp
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (2) :: lrel
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: imcomp
! ..
! .. External Subroutines ..
      EXTERNAL imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imunpk(mb,my)
      ipcomp = imcomp(mx,lrel,my)
      RETURN
    END FUNCTION ipcomp
    SUBROUTINE ipdim(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL imdim, impack, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imunpk(mb,my)
      CALL imdim(mx,my,mx)
      CALL impack(mx,mc)
      RETURN
    END SUBROUTINE ipdim
    SUBROUTINE ipdiv(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL imdiv, impack, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imunpk(mb,my)
      CALL imdiv(mx,my,mx)
      CALL impack(mx,mc)
      RETURN
    END SUBROUTINE ipdiv
    SUBROUTINE ipdivi(ma,ival,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL imdivi, impack, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imdivi(mx,ival,mx)
      CALL impack(mx,mb)
      RETURN
    END SUBROUTINE ipdivi
    SUBROUTINE ipdivr(ma,mb,mc,md)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack), md(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL imdivr, impack, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imunpk(mb,my)
      CALL imdivr(mx,my,mx,my)
      CALL impack(mx,mc)
      CALL impack(my,md)
      RETURN
    END SUBROUTINE ipdivr
    SUBROUTINE ipdvir(ma,ival,mb,irem)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: irem, ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL imdvir, impack, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imdvir(mx,ival,mx,irem)
      CALL impack(mx,mb)
      RETURN
    END SUBROUTINE ipdvir
    SUBROUTINE ipeq(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL imeq, impack, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imeq(mx,my)
      CALL impack(my,mb)
      RETURN
    END SUBROUTINE ipeq
    SUBROUTINE ipfm2i(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmunpk, imfm2i, impack
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL fmunpk(ma,mx)
      CALL imfm2i(mx,mx)
      CALL impack(mx,mb)
      RETURN
    END SUBROUTINE ipfm2i
    SUBROUTINE ipform(form,ma,string)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (*) :: form, string
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL imform, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imform(form,mx,string)
      RETURN
    END SUBROUTINE ipform
    SUBROUTINE ipfprt(form,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (*) :: form
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL imfprt, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imfprt(form,mx)
      RETURN
    END SUBROUTINE ipfprt
    SUBROUTINE ipgcd(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL imgcd, impack, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imunpk(mb,my)
      CALL imgcd(mx,my,mx)
      CALL impack(mx,mc)
      RETURN
    END SUBROUTINE ipgcd
    SUBROUTINE ipi2fm(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack, imi2fm, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imi2fm(mx,mx)
      CALL fmpack(mx,mb)
      RETURN
    END SUBROUTINE ipi2fm
    SUBROUTINE ipi2m(ival,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL imi2m, impack
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imi2m(ival,mx)
      CALL impack(mx,ma)
      RETURN
    END SUBROUTINE ipi2m
    SUBROUTINE ipinp(line,ma,la,lb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: la, lb
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
      CHARACTER (1) :: line(lb)
! ..
! .. External Subroutines ..
      EXTERNAL iminp, impack
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL iminp(line,mx,la,lb)
      CALL impack(mx,ma)
      RETURN
    END SUBROUTINE ipinp
    SUBROUTINE ipm2dp(ma,dval)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL (KIND(0.0D0)) :: dval
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL imm2dp, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imm2dp(mx,dval)
      RETURN
    END SUBROUTINE ipm2dp
    SUBROUTINE ipm2i(ma,ival)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL imm2i, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imm2i(mx,ival)
      RETURN
    END SUBROUTINE ipm2i
    SUBROUTINE ipmax(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL immax, impack, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imunpk(mb,my)
      CALL immax(mx,my,mx)
      CALL impack(mx,mc)
      RETURN
    END SUBROUTINE ipmax
    SUBROUTINE ipmin(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL immin, impack, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imunpk(mb,my)
      CALL immin(mx,my,mx)
      CALL impack(mx,mc)
      RETURN
    END SUBROUTINE ipmin
    SUBROUTINE ipmod(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL immod, impack, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imunpk(mb,my)
      CALL immod(mx,my,mx)
      CALL impack(mx,mc)
      RETURN
    END SUBROUTINE ipmod
    SUBROUTINE ipmpy(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL immpy, impack, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imunpk(mb,my)
      CALL immpy(mx,my,mx)
      CALL impack(mx,mc)
      RETURN
    END SUBROUTINE ipmpy
    SUBROUTINE ipmpyi(ma,ival,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL immpyi, impack, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL immpyi(mx,ival,mx)
      CALL impack(mx,mb)
      RETURN
    END SUBROUTINE ipmpyi
    SUBROUTINE ipmpym(ma,mb,mc,md)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack), md(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL immpym, impack, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imunpk(mb,my)
      CALL imunpk(mc,mz)
      CALL immpym(mx,my,mz,mz)
      CALL impack(mz,md)
      RETURN
    END SUBROUTINE ipmpym
    SUBROUTINE ipout(ma,line,lb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: lb
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
      CHARACTER (1) :: line(lb)
! ..
! .. External Subroutines ..
      EXTERNAL imout, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imout(mx,line,lb)
      RETURN
    END SUBROUTINE ipout
    SUBROUTINE ippmod(ma,mb,mc,md)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack), md(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL impack, impmod, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imunpk(mb,my)
      CALL imunpk(mc,mz)
      CALL impmod(mx,my,mz,mx)
      CALL impack(mx,md)
      RETURN
    END SUBROUTINE ippmod
    SUBROUTINE ipprnt(ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL imprnt, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imprnt(mx)
      RETURN
    END SUBROUTINE ipprnt
    SUBROUTINE ippwr(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL impack, impwr, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imunpk(mb,my)
      CALL impwr(mx,my,mx)
      CALL impack(mx,mc)
      RETURN
    END SUBROUTINE ippwr
    SUBROUTINE ipread(kread,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kread
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL impack, imread
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imread(kread,mx)
      CALL impack(mx,ma)
      RETURN
    END SUBROUTINE ipread
    SUBROUTINE ipsign(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL impack, imsign, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imunpk(mb,my)
      CALL imsign(mx,my,mx)
      CALL impack(mx,mc)
      RETURN
    END SUBROUTINE ipsign
    SUBROUTINE ipsqr(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL impack, imsqr, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imsqr(mx,my)
      CALL impack(my,mb)
      RETURN
    END SUBROUTINE ipsqr
    SUBROUTINE ipst2m(string,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (*) :: string
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL impack, imst2m
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imst2m(string,mx)
      CALL impack(mx,ma)
      RETURN
    END SUBROUTINE ipst2m
    SUBROUTINE ipsub(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack), mb(0:lpack), mc(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL impack, imsub, imunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imunpk(mb,my)
      CALL imsub(mx,my,mx)
      CALL impack(mx,mc)
      RETURN
    END SUBROUTINE ipsub
    SUBROUTINE ipwrit(kwrite,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kwrite
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL imunpk, imwrit
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpck), my(0:lunpck), mz(0:lunpck)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mx, my, mz
! ..

      CALL imunpk(ma,mx)
      CALL imwrit(kwrite,mx)
      RETURN
!             End of the FM package.
    END SUBROUTINE ipwrit
