
!     ZM 1.1                  David M. Smith               5-19-97


!  The ZM routines perform complex floating-point multiple-precision
!  arithmetic.

!  These routines use the FMLIB package (version 1.1) for real
!  floating-point multiple-precision arithmetic.
!  FMLIB 1.0 is Algorithm 693, ACM Transactions on Mathematical
!  Software, Vol. 17, No. 2, June 1991, pages 273-283.

!  This package and FMLIB 1.1 use double precision arithmetic and arrays
!  internally.  This is usually faster at higher precision, and on many
!  machines it is also faster at lower precision.  Both packages are
!  written so that the arithmetic used can easily be changed from double
!  precision to integer, or another available arithmetic type.  See the
!  EFFICIENCY discussion in FMLIB 1.1 for details.



!  1. INITIALIZING THE PACKAGE

!  Before calling any routine in the package, several variables in the
!  common blocks /FMUSER/, /FM/, /FMSAVE/, /FMBUFF/, and /ZMUSER/ must
!  be initialized.  These common blocks contain information that is
!  saved between calls, so they should be declared in the main program.

!  Subroutine ZMSET initializes these variables to default values and
!  defines all machine-dependent values in the package.  After calling
!  ZMSET once at the start of a program, the user may sometimes want to
!  reset some of the variables in common blocks /FMUSER/ or /ZMUSER/.


!  2.  REPRESENTATION OF ZM NUMBERS

!  The format for complex FM numbers (called ZM numbers below) is very
!  similar to that for real FM numbers in FMLIB.  Each ZM array holds
!  two FM numbers to represent the real and imaginary parts of a complex
!  number.  Each ZM array is twice as long as a corresponding FM array,
!  with the imaginary part starting at the midpoint of the array.  As
!  with FM, there are packed and unpacked formats for the numbers.


!  3. INPUT/OUTPUT ROUTINES

!  All versions of the input routines perform free-format conversion
!  from characters to ZM numbers.

!  a. Conversion to or from a character array

!     ZMINP converts from a character*1 array to an ZM number.

!     ZMOUT converts an ZM number to base 10 and formats it for output
!           as an array of type character*1.  The output is left
!           justified in the array, and the format is defined by
!           variables in common, so that a separate format definition
!           does not have to be provided for each output call.

!     For the output format of ZM numbers, JFORM1 and JFORM2 determine
!     the format for the individual parts of a complex number as
!     described in the FMLIB documentation.

!     JFORMZ (in /ZMUSER/) determines the combined output format of the
!     real and imaginary parts.

!     JFORMZ = 1  normal setting :       1.23 - 4.56 i
!            = 2  use capital I  :       1.23 - 4.56 I
!            = 3  parenthesis format   ( 1.23 , -4.56 )

!     JPRNTZ (in /ZMUSER/) controls whether to print real
!     and imaginary parts on one line whenever possible.

!     JPRNTZ = 1  print both parts as a single string :
!                     1.23456789M+321 - 9.87654321M-123 i
!            = 2  print on separate lines without the 'i' :
!                     1.23456789M+321
!                    -9.87654321M-123

!  b. Conversion to or from a character string

!     ZMST2M converts from a character string to an ZM number.

!     ZMFORM converts an ZM number to a character string according to
!            a format provided in each call.  The format descriptions
!            are more like that of a Fortran FORMAT statement, and
!            integer or fixed-point output is right justified.

!  c. Direct read or write

!     ZMPRNT uses ZMOUT to print one ZM number.

!     ZMFPRT uses ZMFORM to print one ZM number.

!     ZMWRIT writes ZM numbers for later input using ZMREAD.

!     ZMREAD reads ZM numbers written by ZMWRIT.

!  For further description of these routines, see section 5 below.


!  4. ARRAY DIMENSIONS

!  The parameters LPACKZ and LUNPKZ define the size of the packed and
!  unpacked ZM arrays.  The real part starts at the beginning of the
!  array, and the imaginary part starts at word KPTIMP for packed format
!  or at word KPTIMU for unpacked format.


!  5.  LIST OF ROUTINES

!  These are the routines in ZMLIB that are designed to be called by
!  the user.  All are subroutines, and in each case the version of the
!  routine to handle packed ZM numbers has the same name, with 'ZM'
!  replaced by 'ZP'.

!  MA, MB, MC refer to ZM format complex numbers.
!  MAFM, MBFM, MCFM refer to FM format real numbers.
!  INTEG is a Fortran INTEGER variable.
!  ZVAL is a Fortran COMPLEX variable.

!  In each case it is permissible to use the same array more than
!  once in the calling sequence.  The statement
!  MA = MA*MA   may be written  CALL ZMMPY(MA,MA,MA).

!  ZMABS(MA,MBFM)       MBFM = ABS(MA)    Result is real.

!  ZMACOS(MA,MB)        MB = ACOS(MA)

!  ZMADD(MA,MB,MC)      MC = MA + MB

!  ZMADDI(MA,INTEG)     MA = MA + INTEG  Increment an ZM number by a one
!                                        word integer.  Note this call
!                                        does not have an "MB" result
!                                        like ZMDIVI and ZMMPYI.

!  ZMARG(MA,MBFM)       MBFM = Argument(MA)    Result is real.

!  ZMASIN(MA,MB)        MB = ASIN(MA)

!  ZMATAN(MA,MB)        MB = ATAN(MA)

!  ZMCHSH(MA,MB,MC)     MB = COSH(MA),  MC = SINH(MA).
!                            Faster than 2 calls.

!  ZMCMPX(MAFM,MBFM,MC) MC = CMPLX(MAFM,MBFM)

!  ZMCONJ(MA,MB)        MB = CONJG(MA)

!  ZMCOS(MA,MB)         MB = COS(MA)

!  ZMCOSH(MA,MB)        MB = COSH(MA)

!  ZMCSSN(MA,MB,MC)     MB = COS(MA),  MC = SIN(MA).
!                            Faster than 2 calls.

!  ZMDIV(MA,MB,MC)      MC = MA / MB

!  ZMDIVI(MA,INTEG,MB)  MB = MA / INTEG

!  ZMEQ(MA,MB)          MB = MA

!  ZMEQU(MA,MB,NDA,NDB) MB = MA    Version for changing precision.
!                                  (NDA and NDB are as in FMEQU)

!  ZMEXP(MA,MB)         MB = EXP(MA)

!  ZMFORM(FORM1,FORM2,MA,STRING)   STRING = MA
!                       MA is converted to a character string using
!                       format FORM1 for the real part and FORM2 for
!                       the imaginary part.  The  result is returned
!                       in STRING.  FORM1 and FORM2 can represent I,
!                       F, E, or 1PE formats.  Example:
!                          CALL ZMFORM('F20.10','F15.10',MA,STRING)

!  ZMFPRT(FORM1,FORM2,MA)    Print MA on unit KW using
!                            formats FORM1 and FORM2.

!  ZMI2M(INTEG,MA)           MA = CMPLX(INTEG,0)

!  ZM2I2M(INTEG1,INTEG2,MA)  MA = CMPLX(INTEG1,INTEG2)

!  ZMIMAG(MA,MBFM)           MBFM = IMAG(MA)    Imaginary part.

!  ZMINP(LINE,MA,LA,LB)      MA = LINE   Input conversion.
!                                 Convert LINE(LA) through LINE(LB)
!                                 from characters to ZM.  LINE is a
!                                 character array of length at least LB.

!  ZMINT(MA,MB)         MB = INT(MA)        Integer part of both Real
!                                           and Imaginary parts of MA.

!  ZMIPWR(MA,INTEG,MB)  MB = MA ** INTEG    Integer power function.

!  ZMLG10(MA,MB)        MB = LOG10(MA)

!  ZMLN(MA,MB)          MB = LOG(MA)

!  ZMM2I(MA,INTEG)      INTEG = INT(REAL(MA))

!  ZMM2Z(MA,ZVAL)       ZVAL = MA

!  ZMMPY(MA,MB,MC)      MC = MA * MB

!  ZMMPYI(MA,INTEG,MB)  MB = MA * INTEG

!  ZMNINT(MA,MB)        MB = NINT(MA)   Nearest integer of both Real
!                                       and Imaginary.

!  ZMOUT(MA,LINE,LB,LAST1,LAST2)        LINE = MA
!                       Convert from FM to character.
!                       LINE  is the returned character array.
!                       LB    is the dimensioned size of LINE.
!                       LAST1 is returned as the position in LINE of
!                             the last character of REAL(MA).
!                       LAST2 is returned as the position in LINE
!                             of the last character of AIMAG(MA).

!  ZMPRNT(MA)           Print MA on unit KW using current format.

!  ZMPWR(MA,MB,MC)      MC = MA ** MB

!  ZMREAD(KREAD,MA)     MA   is returned after reading one (possibly
!                            multi-line) ZM number on unit KREAD.  This
!                            routine reads numbers written by ZMWRIT.

!  ZMREAL(MA,MBFM)      MBFM = REAL(MA)    Real part.

!  ZMRPWR(MA,IVAL,JVAL,MB)     MB = MA ** (IVAL/JVAL)

!  ZMSET(NPREC)         Initialize ZM package.  Set precision to the
!                       equivalent of at least NPREC base 10 digits.

!  ZMSIN(MA,MB)         MB = SIN(MA)

!  ZMSINH(MA,MB)        MB = SINH(MA)

!  ZMSQR(MA,MB)         MB = MA*MA    Faster than ZMMPY.

!  ZMSQRT(MA,MB)        MB = SQRT(MA)

!  ZMST2M(STRING,MA)    MA = STRING
!                            Convert from character string to ZM.
!                            Often more convenient than ZMINP, which
!                            converts an array of CHARACTER*1 values.
!                            Example: CALL ZMST2M('123.4+5.67i',MA).

!  ZMSUB(MA,MB,MC)      MC = MA - MB

!  ZMTAN(MA,MB)         MB = TAN(MA)

!  ZMTANH(MA,MB)        MB = TANH(MA)

!  ZMWRIT(KWRITE,MA)    Write MA on unit KWRITE.  Multi-line numbers
!                       are formatted for automatic reading with ZMREAD.

!  ZMZ2M(ZVAL,MA)       MA = ZVAL



    SUBROUTINE zmset(nprec)

!  Initialize common and set precision to at least NPREC significant
!  digits before using ZM arithmetic.

      IMPLICIT NONE

!  Here are the common blocks used for complex arithmetic.

!  /FMUSER/, /FM/, /FMBUFF/, /FMSAVE/, and /ZMUSER/ should also be
!  declared in the main program.

! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: nprec
! ..
! .. External Subroutines ..
      EXTERNAL fmset
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, dlogtw, dpeps, &
        dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli, mbspi, mexpab, &
        mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, jformz, jprntz, kaccsw, &
        kdebug, keswch, kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, &
        ncall, ndg2mx, ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, &
        ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa), mx(0:lunpkz), my(0:lunpkz), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff), cmbufz(lmbufz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
      COMMON /zmbuff/cmbufz
      COMMON /zmpck/mx, my
      COMMON /zmuser/jformz, jprntz
! ..
!             Set JFORMZ to ' 1.23 + 4.56 i ' format.

      jformz = 1

!             Set JPRNTZ to print real and imaginary parts on one
!             line whenever possible.

      jprntz = 1

!             Use FMSET to initialize the other common blocks.

      CALL fmset(nprec)

      RETURN
    END SUBROUTINE zmset
    SUBROUTINE zmabs(ma,mbfm)

!  MBFM = ABS(MA)

!  Complex absolute value.  The result is a real FM number.

      IMPLICIT NONE

!             Scratch array usage during ZMABS:   M01 - M02, MZ01

! .. Intrinsic Functions ..
      INTRINSIC int, max, min
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mbfm(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccmb, maiz, marz, mxexp1, mxsave
      INTEGER :: kasave, kovun, kreslt, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmabs, fmadd, fmeq, fmi2m, fmsqr, fmsqrt, zmentr, zmeq2, zmexi2
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
        m05(0:lunpck), m06(0:lunpck), mwa(lmwa), mz01(0:lunpkz), &
        mz02(0:lunpkz), mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      CALL zmentr('ZMABS ',ma,ma,1,mz01,kreslt,ndsave,mxsave,kasave,kovun)
      IF (kreslt/=0) THEN
        CALL fmeq(mz01,mbfm)
        RETURN
      END IF
      marz = ma(0)
      maiz = ma(kptimu)
      kaccsw = 0
      CALL zmeq2(ma,ma,ndsave,ndig,1)

!             Check for special cases.

      mxexp1 = int(mxexp2/2.01D0)
      IF (ma(2)==0) THEN
        CALL fmabs(ma(kptimu),mbfm)
        GO TO 10
      ELSE IF (ma(kptimu+2)==0) THEN
        CALL fmabs(ma,mbfm)
        GO TO 10
      ELSE IF (ma(1)==mexpov .OR. ma(kptimu+1)==mexpov) THEN
        CALL fmi2m(1,mbfm)
        mbfm(1) = max(ma(1),ma(kptimu+1))
        GO TO 10
      ELSE IF (ma(1)==mexpun) THEN
        IF (ma(kptimu+1)>-mxexp1+ndig+1) THEN
          CALL fmabs(ma(kptimu),mbfm)
          GO TO 10
        END IF
      ELSE IF (ma(kptimu+1)==mexpun) THEN
        IF (ma(1)>-mxexp1+ndig+1) THEN
          CALL fmabs(ma,mbfm)
          GO TO 10
        END IF
      ELSE IF (ma(1)/=munkno .AND. ma(kptimu+1)/=munkno) THEN
        IF (ma(1)>ma(kptimu+1)+ndig+1) THEN
          CALL fmabs(ma,mbfm)
          GO TO 10
        ELSE IF (ma(kptimu+1)>ma(1)+ndig+1) THEN
          CALL fmabs(ma(kptimu),mbfm)
          GO TO 10
        END IF
      END IF

      CALL fmsqr(ma,m01)
      CALL fmsqr(ma(kptimu),m02)
      CALL fmadd(m01,m02,mbfm)
      CALL fmsqrt(mbfm,mbfm)

10    maccmb = mbfm(0)
      ma(0) = marz
      ma(kptimu) = maiz
      mbfm(0) = min(maccmb,marz,maiz)
      CALL zmexi2(mbfm,mbfm,ndsave,mxsave,kasave,kovun,1)
      RETURN
    END SUBROUTINE zmabs
    SUBROUTINE zmacos(ma,mb)

!  MB = ACOS(MA).

      IMPLICIT NONE

!             Scratch array usage during ZMACOS:  M01 - M06, MZ01 - MZ03

! .. Intrinsic Functions ..
      INTRINSIC int, min
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccmb, maiz, marz, mxsave
      INTEGER :: j, kasave, kovun, kreslt, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmacos, fmadd, fmdivi, fmi2m, fmpi, fmsqr, fmsub, zmadd, &
        zmentr, zmeq2, zmexit, zmi2m, zmln, zmmpy, zmntr, zmrslt, zmsqrt, &
        zmsub, zmwarn
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
        mpisav(0:lunpck), mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      CALL zmentr('ZMACOS',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
      IF (kreslt/=0) RETURN
      marz = ma(0)
      maiz = ma(kptimu)
      kaccsw = 0

      CALL zmeq2(ma,ma,ndsave,ndig,1)

!             Check for special cases.

      IF (ma(2)==0 .AND. ma(kptimu+2)==0) THEN
        CALL fmpi(mz01)
        CALL fmdivi(mz01,2,mz01)
        CALL fmi2m(0,mz01(kptimu))
        GO TO 60
      ELSE IF (ma(kptimu+2)==0) THEN
        CALL fmacos(ma,mz01)
        IF (kflag==0) THEN
          CALL fmi2m(0,mz01(kptimu))
          GO TO 60
        END IF
      END IF
      IF ((ma(2)==0 .OR. ma(1)*2<=-ndig) .AND. (ma(kptimu+ &
          2)==0 .OR. ma(kptimu+1)*2<=-ndig)) THEN
        CALL fmpi(mz01)
        CALL fmdivi(mz01,2,mz01)
        CALL fmi2m(0,mz01(kptimu))
        CALL zmsub(mz01,ma,mz01)
        GO TO 60
      END IF

      CALL zmi2m(1,mz03)
      CALL zmsub(mz03,ma,mz02)
      CALL zmadd(mz03,ma,mz03)
      CALL zmmpy(mz02,mz03,mz02)
      CALL zmsqrt(mz02,mz02)
      DO 10 j = 0, ndig + 1
        mz03(j) = mz02(kptimu+j)
        mz03(kptimu+j) = mz02(j)
10    CONTINUE
      IF (mz03(1)/=munkno) mz03(2) = -mz03(2)

      IF ((ma(2)/=0 .AND. mz03(1)==ma(1) .AND. mz03(2)==ma(2)) .OR. (ma( &
          kptimu+2)/=0 .AND. mz03(kptimu+1)==ma(kptimu+1) .AND. mz03(kptimu+ &
          2)==ma(kptimu+2))) THEN
        CALL zmadd(ma,mz03,mz03)

        CALL fmsqr(mz03,m04)
        CALL fmsqr(mz03(kptimu),m05)
        CALL fmadd(m04,m05,m06)
        CALL fmi2m(1,m03)
        CALL fmsub(m06,m03,m03)
        IF (m03(1)<0) THEN
          ndig = ndig - int(m03(1))
          IF (ndig>ndg2mx) THEN
            namest(ncall) = 'ZMACOS'
            kflag = -9
            CALL zmwarn
            kreslt = 12
            ndig = ndsave
            CALL zmrslt(mb,kreslt)
            IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
            ncall = ncall - 1
            mxexp = mxsave
            kaccsw = kasave
            RETURN
          END IF
          CALL zmeq2(ma,ma,ndsave,ndig,1)
          CALL zmi2m(1,mz03)
          CALL zmsub(mz03,ma,mz02)
          CALL zmadd(mz03,ma,mz03)
          CALL zmmpy(mz02,mz03,mz02)
          CALL zmsqrt(mz02,mz02)
          DO 20 j = 0, ndig + 1
            mz03(j) = mz02(kptimu+j)
            mz03(kptimu+j) = mz02(j)
20        CONTINUE
          IF (mz03(1)/=munkno) mz03(2) = -mz03(2)
          CALL zmadd(ma,mz03,mz03)
        END IF

        CALL zmln(mz03,mz03)
        DO 30 j = 0, ndig + 1
          mz01(j) = mz03(kptimu+j)
          mz01(kptimu+j) = mz03(j)
30      CONTINUE
        IF (mz01(kptimu+1)/=munkno) mz01(kptimu+2) = -mz01(kptimu+2)
      ELSE
        CALL zmsub(ma,mz03,mz03)

        CALL fmsqr(mz03,m04)
        CALL fmsqr(mz03(kptimu),m05)
        CALL fmadd(m04,m05,m06)
        CALL fmi2m(1,m03)
        CALL fmsub(m06,m03,m03)
        IF (m03(1)<0) THEN
          ndig = ndig - int(m03(1))
          IF (ndig>ndg2mx) THEN
            namest(ncall) = 'ZMACOS'
            kflag = -9
            CALL zmwarn
            kreslt = 12
            ndig = ndsave
            CALL zmrslt(mb,kreslt)
            IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
            ncall = ncall - 1
            mxexp = mxsave
            kaccsw = kasave
            RETURN
          END IF
          CALL zmeq2(ma,ma,ndsave,ndig,1)
          CALL zmi2m(1,mz03)
          CALL zmsub(mz03,ma,mz02)
          CALL zmadd(mz03,ma,mz03)
          CALL zmmpy(mz02,mz03,mz02)
          CALL zmsqrt(mz02,mz02)
          DO 40 j = 0, ndig + 1
            mz03(j) = mz02(kptimu+j)
            mz03(kptimu+j) = mz02(j)
40        CONTINUE
          IF (mz03(1)/=munkno) mz03(2) = -mz03(2)
          CALL zmsub(ma,mz03,mz03)
        END IF

        CALL zmln(mz03,mz03)
        DO 50 j = 0, ndig + 1
          mz01(j) = mz03(kptimu+j)
          mz01(kptimu+j) = mz03(j)
50      CONTINUE
        IF (mz01(1)/=munkno) mz01(2) = -mz01(2)
      END IF

60    maccmb = mz01(0)
      ma(0) = marz
      mz01(0) = min(maccmb,marz,maiz)
      maccmb = mz01(kptimu)
      ma(kptimu) = maiz
      mz01(kptimu) = min(maccmb,marz,maiz)
      CALL zmexit(mz01,mb,ndsave,mxsave,kasave,kovun,0)
      RETURN
    END SUBROUTINE zmacos
    SUBROUTINE zmadd(ma,mb,mc)

!  MC = MA + MB

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz), mc(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: mai, mar, mbi, mbr, mxsave
      INTEGER :: kasave, kf1, kovun, kreslt, kwrnsv, ndsave, ntrsav
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, zmentr, zmntr, zmwarn
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
      IF (abs(ma(1))>mexpab .OR. abs(ma(kptimu+1))>mexpab .OR. abs(mb(1))> &
          mexpab .OR. abs(mb(kptimu+1))>mexpab .OR. kdebug>=1) THEN
        CALL zmentr('ZMADD ',ma,mb,2,mc,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
        ndig = ndsave
        mxexp = mxsave
        kaccsw = kasave
        ntrsav = ntrace
      ELSE
        ncall = ncall + 1
        ntrsav = ntrace
        IF (ntrace/=0) THEN
          namest(ncall) = 'ZMADD '
          CALL zmntr(2,ma,mb,2)
          ntrace = 0
        END IF
        kovun = 0
      END IF

!             Force FMADD to use more guard digits for user calls.

      ncall = ncall - 1

      kwrnsv = kwarn
      kwarn = 0
      mar = ma(1)
      IF (ma(2)==0) mar = mexpun - 1
      mai = ma(kptimu+1)
      IF (ma(kptimu+2)==0) mai = mexpun - 1
      mbr = mb(1)
      IF (mb(2)==0) mbr = mexpun - 1
      mbi = mb(kptimu+1)
      IF (mb(kptimu+2)==0) mbi = mexpun - 1

      CALL fmadd(ma,mb,mc)
      kf1 = kflag
      CALL fmadd(ma(kptimu),mb(kptimu),mc(kptimu))

      ncall = ncall + 1
      IF (ntrsav/=0) THEN
        ntrace = ntrsav
        namest(ncall) = 'ZMADD '
      END IF
      kwarn = kwrnsv
      IF (kflag==1) kflag = kf1
      IF (kflag==1) THEN
        kflag = 0
        IF (mar<=mbr .AND. mai<=mbi) kflag = 1
        IF (mar>=mbr .AND. mai>=mbi) kflag = 1
      END IF

      IF (mc(1)==munkno .OR. mc(kptimu+1)==munkno) THEN
        kflag = -4
      ELSE IF (mc(1)==mexpov .OR. mc(kptimu+1)==mexpov) THEN
        kflag = -5
      ELSE IF (mc(1)==mexpun .OR. mc(kptimu+1)==mexpun) THEN
        kflag = -6
      END IF
      IF ((mc(1)==munkno) .OR. (mc(kptimu+1)==munkno) .OR. (mc(1)==mexpun &
          .AND. kovun==0) .OR. (mc(kptimu+1)==mexpun .AND. kovun==0) .OR. (mc( &
          1)==mexpov .AND. kovun==0) .OR. (mc(kptimu+ &
          1)==mexpov .AND. kovun==0)) THEN
        namest(ncall) = 'ZMADD '
        CALL zmwarn
      END IF
      IF (ntrace/=0) THEN
        CALL zmntr(1,mc,mc,1)
      END IF
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmadd
    SUBROUTINE zmaddi(ma,integ)

!  MA = MA + INTEG        Increment by one-word (real) integer.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: integ
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: mxsave
      INTEGER :: kasave, kovun, kreslt, kwrnsv, ndsave, ntrsav
! ..
! .. External Subroutines ..
      EXTERNAL fmaddi, fmntri, zmentr, zmntr, zmwarn
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
      IF (abs(ma(1))>mexpab .OR. abs(ma(kptimu+1))>mexpab .OR. kdebug>=1) THEN
        CALL zmentr('ZMADDI',ma,ma,1,ma,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
        ndig = ndsave
        mxexp = mxsave
        kaccsw = kasave
        ntrsav = ntrace
      ELSE
        ncall = ncall + 1
        IF (ntrace/=0) THEN
          namest(ncall) = 'ZMADDI'
          CALL zmntr(2,ma,ma,1)
          CALL fmntri(2,integ,0)
        END IF
        kovun = 0
      END IF

!             Force FMADDI to use more guard digits for user calls.

      ncall = ncall - 1
      ntrsav = ntrace
      ntrace = 0
      kwrnsv = kwarn
      kwarn = 0

      CALL fmaddi(ma,integ)

      ntrace = ntrsav
      kwarn = kwrnsv
      ncall = ncall + 1
      IF (ntrace/=0) namest(ncall) = 'ZMADDI'
      IF (ma(1)==munkno .OR. ma(kptimu+1)==munkno) THEN
        kflag = -4
      ELSE IF (ma(1)==mexpov .OR. ma(kptimu+1)==mexpov) THEN
        kflag = -5
      ELSE IF (ma(1)==mexpun .OR. ma(kptimu+1)==mexpun) THEN
        kflag = -6
      END IF
      IF ((ma(1)==munkno) .OR. (ma(kptimu+1)==munkno) .OR. (ma(1)==mexpun &
          .AND. kovun==0) .OR. (ma(kptimu+1)==mexpun .AND. kovun==0) .OR. (ma( &
          1)==mexpov .AND. kovun==0) .OR. (ma(kptimu+ &
          1)==mexpov .AND. kovun==0)) THEN
        namest(ncall) = 'ZMADDI'
        CALL zmwarn
      END IF
      IF (ntrace/=0) CALL zmntr(1,ma,ma,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmaddi
    SUBROUTINE zmarg(ma,mbfm)

!  MBFM = ARG(MA)

!  Complex argument.  The result is a real FM number.

      IMPLICIT NONE

!             Scratch array usage during ZMARG:   M01 - M06, MZ01

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mbfm(0:lunpck)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: mxsave
      INTEGER :: kasave, kovun, kreslt, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmatn2, fmeq, zmentr, zmeq2, zmexi2
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
        m05(0:lunpck), m06(0:lunpck), mwa(lmwa), mz01(0:lunpkz), &
        mz02(0:lunpkz), mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      CALL zmentr('ZMARG ',ma,ma,1,mz01,kreslt,ndsave,mxsave,kasave,kovun)
      IF (kreslt/=0) THEN
        CALL fmeq(mz01,mbfm)
        RETURN
      END IF
      kaccsw = 0
      CALL zmeq2(ma,ma,ndsave,ndig,1)

      CALL fmatn2(ma(kptimu),ma,mbfm)

      CALL zmexi2(mbfm,mbfm,ndsave,mxsave,kasave,kovun,1)
      RETURN
    END SUBROUTINE zmarg
    SUBROUTINE zmasin(ma,mb)

!  MB = ASIN(MA).

      IMPLICIT NONE

!             Scratch array usage during ZMASIN:  M01 - M06, MZ01 - MZ03

! .. Intrinsic Functions ..
      INTRINSIC int, min
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccmb, maiz, marz, mxsave
      INTEGER :: j, kasave, kovun, kreslt, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmasin, fmi2m, fmsqr, fmsub, zmadd, zmentr, zmeq, zmeq2, &
        zmexit, zmi2m, zmln, zmmpy, zmntr, zmrslt, zmsqrt, zmsub, zmwarn
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
        mpisav(0:lunpck), mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      CALL zmentr('ZMASIN',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
      IF (kreslt/=0) RETURN
      marz = ma(0)
      maiz = ma(kptimu)
      kaccsw = 0

      CALL zmeq2(ma,ma,ndsave,ndig,1)

!             Check for special cases.

      IF (ma(2)==0 .AND. ma(kptimu+2)==0) THEN
        CALL zmi2m(0,mz01)
        GO TO 60
      ELSE IF ((ma(2)==0 .OR. ma(1)*2<=-ndig) .AND. (ma(kptimu+ &
          2)==0 .OR. ma(kptimu+1)*2<=-ndig)) THEN
        CALL zmeq(ma,mz01)
        GO TO 60
      ELSE IF (ma(kptimu+2)==0) THEN
        CALL fmasin(ma,mz01)
        IF (kflag==0) THEN
          CALL fmi2m(0,mz01(kptimu))
          GO TO 60
        END IF
      END IF

      CALL zmi2m(1,mz03)
      CALL zmsub(mz03,ma,mz02)
      CALL zmadd(mz03,ma,mz03)
      CALL zmmpy(mz02,mz03,mz02)
      CALL zmsqrt(mz02,mz02)
      DO 10 j = 0, ndig + 1
        mz03(j) = ma(kptimu+j)
        mz03(kptimu+j) = ma(j)
10    CONTINUE
      IF (mz03(1)/=munkno) mz03(2) = -mz03(2)

      IF ((mz02(2)/=0 .AND. mz03(1)==mz02(1) .AND. mz03(2)==mz02( &
          2)) .OR. (mz02(kptimu+2)/=0 .AND. mz03(kptimu+1)==mz02(kptimu+ &
          1) .AND. mz03(kptimu+2)==mz02(kptimu+2))) THEN
        CALL zmadd(mz02,mz03,mz03)

        CALL fmsqr(mz03,m04)
        CALL fmsqr(mz03(kptimu),m05)
        CALL fmadd(m04,m05,m06)
        CALL fmi2m(1,m03)
        CALL fmsub(m06,m03,m03)
        IF (m03(1)<0) THEN
          ndig = ndig - int(m03(1))
          IF (ndig>ndg2mx) THEN
            namest(ncall) = 'ZMASIN'
            kflag = -9
            CALL zmwarn
            kreslt = 12
            ndig = ndsave
            CALL zmrslt(mb,kreslt)
            IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
            ncall = ncall - 1
            mxexp = mxsave
            kaccsw = kasave
            RETURN
          END IF
          CALL zmeq2(ma,ma,ndsave,ndig,1)
          CALL zmi2m(1,mz03)
          CALL zmsub(mz03,ma,mz02)
          CALL zmadd(mz03,ma,mz03)
          CALL zmmpy(mz02,mz03,mz02)
          CALL zmsqrt(mz02,mz02)
          DO 20 j = 0, ndig + 1
            mz03(j) = ma(kptimu+j)
            mz03(kptimu+j) = ma(j)
20        CONTINUE
          IF (mz03(1)/=munkno) mz03(2) = -mz03(2)
          CALL zmadd(mz02,mz03,mz03)
        END IF

        CALL zmln(mz03,mz03)
        DO 30 j = 0, ndig + 1
          mz01(j) = mz03(kptimu+j)
          mz01(kptimu+j) = mz03(j)
30      CONTINUE
        IF (mz01(kptimu+1)/=munkno) mz01(kptimu+2) = -mz01(kptimu+2)
      ELSE
        CALL zmsub(mz02,mz03,mz03)

        CALL fmsqr(mz03,m04)
        CALL fmsqr(mz03(kptimu),m05)
        CALL fmadd(m04,m05,m06)
        CALL fmi2m(1,m03)
        CALL fmsub(m06,m03,m03)
        IF (m03(1)<0) THEN
          ndig = ndig - int(m03(1))
          IF (ndig>ndg2mx) THEN
            namest(ncall) = 'ZMASIN'
            kflag = -9
            CALL zmwarn
            kreslt = 12
            ndig = ndsave
            CALL zmrslt(mb,kreslt)
            IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
            ncall = ncall - 1
            mxexp = mxsave
            kaccsw = kasave
            RETURN
          END IF
          CALL zmeq2(ma,ma,ndsave,ndig,1)
          CALL zmi2m(1,mz03)
          CALL zmsub(mz03,ma,mz02)
          CALL zmadd(mz03,ma,mz03)
          CALL zmmpy(mz02,mz03,mz02)
          CALL zmsqrt(mz02,mz02)
          DO 40 j = 0, ndig + 1
            mz03(j) = ma(kptimu+j)
            mz03(kptimu+j) = ma(j)
40        CONTINUE
          IF (mz03(1)/=munkno) mz03(2) = -mz03(2)
          CALL zmsub(mz02,mz03,mz03)
        END IF
        CALL zmln(mz03,mz03)
        DO 50 j = 0, ndig + 1
          mz01(j) = mz03(kptimu+j)
          mz01(kptimu+j) = mz03(j)
50      CONTINUE
        IF (mz01(1)/=munkno) mz01(2) = -mz01(2)
      END IF

60    maccmb = mz01(0)
      ma(0) = marz
      mz01(0) = min(maccmb,marz,maiz)
      maccmb = mz01(kptimu)
      ma(kptimu) = maiz
      mz01(kptimu) = min(maccmb,marz,maiz)
      CALL zmexit(mz01,mb,ndsave,mxsave,kasave,kovun,0)
      RETURN
    END SUBROUTINE zmasin
    SUBROUTINE zmatan(ma,mb)

!  MB = ATAN(MA).

      IMPLICIT NONE

!             Scratch array usage during ZMATAN:  M01 - M06, MZ01 - MZ04

! .. Intrinsic Functions ..
      INTRINSIC int, min
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccmb, maiz, marz, mxsave
      REAL :: x
      INTEGER :: j, jterm, kasave, kovun, kreslt, ndsave
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL fmabs, fmadd, fmatan, fmdivi, fmeq, fmi2m, fmpi, fmsp2m, fmsqr, &
        fmsub, zm2i2m, zmadd, zmdiv, zmdivi, zmentr, zmeq, zmeq2, zmexit, &
        zmi2m, zmln, zmmpy, zmntr, zmrslt, zmsqr, zmsub, zmwarn
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
        mpisav(0:lunpck), mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      CALL zmentr('ZMATAN',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
      IF (kreslt/=0) RETURN
      marz = ma(0)
      maiz = ma(kptimu)
      kaccsw = 0

      CALL zmeq2(ma,ma,ndsave,ndig,1)

!             Check for special cases.

      IF (ma(2)==0 .AND. ma(kptimu+2)==0) THEN
        CALL zmi2m(0,mz04)
        GO TO 30
      ELSE IF ((ma(2)==0 .OR. ma(1)*2<=-ndig) .AND. (ma(kptimu+ &
          2)==0 .OR. ma(kptimu+1)*2<=-ndig)) THEN
        CALL zmeq(ma,mz04)
        GO TO 30
      ELSE IF (ma(kptimu+2)==0) THEN
        CALL fmatan(ma,mz04)
        IF (kflag==0) THEN
          CALL fmi2m(0,mz04(kptimu))
          GO TO 30
        END IF
      END IF

      x = 1.0E+5
      CALL fmsp2m(x,m02)
      CALL fmabs(ma,m03)
      CALL fmabs(ma(kptimu),m04)
      CALL fmadd(m03,m04,m04)

      IF (fmcomp(m04,'GE',m02)) THEN
        CALL zmi2m(0,mz04)
        CALL fmpi(mz04)
        CALL fmdivi(mz04,2,mz04)
        IF (ma(2)<0) mz04(2) = -mz04(2)
        CALL zmi2m(1,mz02)
        CALL zmdiv(mz02,ma,mz02)
        CALL zmeq(mz02,mz03)
        CALL zmsub(mz04,mz02,mz04)
        IF (ma(1)>ndig .OR. ma(kptimu+1)>ndig) GO TO 30
        CALL zmsqr(mz02,mz02)
        jterm = 1
10      CALL zmmpy(mz03,mz02,mz03)
        jterm = jterm + 2
        CALL fmeq(mz03,m05)
        CALL fmeq(mz03(kptimu),m06)
        CALL zmdivi(mz03,jterm,mz03)
        CALL zmadd(mz04,mz03,mz04)
        IF (kflag/=0) GO TO 30
        CALL fmeq(m05,mz03)
        CALL fmeq(m06,mz03(kptimu))
        CALL zmmpy(mz03,mz02,mz03)
        jterm = jterm + 2
        CALL fmeq(mz03,m05)
        CALL fmeq(mz03(kptimu),m06)
        CALL zmdivi(mz03,jterm,mz03)
        CALL zmsub(mz04,mz03,mz04)
        IF (kflag/=0) GO TO 30
        CALL fmeq(m05,mz03)
        CALL fmeq(m06,mz03(kptimu))
        GO TO 10
      ELSE
        CALL zm2i2m(0,1,mz02)
        CALL zmsub(mz02,ma,mz03)
        CALL zmadd(mz02,ma,mz02)
        CALL zmdiv(mz02,mz03,mz03)

        CALL fmsqr(mz03,m04)
        CALL fmsqr(mz03(kptimu),m05)
        CALL fmadd(m04,m05,m06)
        CALL fmi2m(1,m03)
        CALL fmsub(m06,m03,m03)
        IF (m03(1)<0) THEN
          ndig = ndig - int(m03(1))
          IF (ndig>ndg2mx) THEN
            namest(ncall) = 'ZMATAN'
            kflag = -9
            CALL zmwarn
            kreslt = 12
            ndig = ndsave
            CALL zmrslt(mb,kreslt)
            IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
            ncall = ncall - 1
            mxexp = mxsave
            kaccsw = kasave
            RETURN
          END IF
          CALL zmeq2(ma,ma,ndsave,ndig,1)
          CALL zm2i2m(0,1,mz02)
          CALL zmsub(mz02,ma,mz03)
          CALL zmadd(mz02,ma,mz02)
          CALL zmdiv(mz02,mz03,mz03)
        END IF

        CALL zmln(mz03,mz03)
        CALL zmdivi(mz03,2,mz03)
        DO 20 j = 0, ndig + 1
          mz04(j) = mz03(kptimu+j)
          mz04(kptimu+j) = mz03(j)
20      CONTINUE
        IF (mz04(1)/=munkno) mz04(2) = -mz04(2)
      END IF

30    maccmb = mz04(0)
      ma(0) = marz
      mz04(0) = min(maccmb,marz,maiz)
      maccmb = mz04(kptimu)
      ma(kptimu) = maiz
      mz04(kptimu) = min(maccmb,marz,maiz)
      CALL zmexit(mz04,mb,ndsave,mxsave,kasave,kovun,0)
      RETURN
    END SUBROUTINE zmatan
    SUBROUTINE zmchsh(ma,mb,mc)

!  MB = COSH(MA),    MC = SINH(MA).

!  If both the hyperbolic sine and cosine are needed, this routine
!  is faster than calling both ZMCOS and ZMSIN.

!  MB and MC must be distinct arrays.

      IMPLICIT NONE

!             Scratch array usage during ZMCHSH:  M01 - M06, MZ01 - MZ04

! .. Intrinsic Functions ..
      INTRINSIC abs, min
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz), mc(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccmb, maiz, marz, mxsave
      INTEGER :: kasave, kovun, kreslt, krsave, ncsave, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmchsh, fmcssn, fmi2m, fmmpy, zmentr, zmeq, zmeq2, zmexit, &
        zmi2m, zmntr, zmntrj, zmprnt
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
        mpisav(0:lunpck), mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      ncsave = ncall
      CALL zmentr('ZMCHSH',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
      ncall = ncsave + 1
      IF (kreslt/=0) THEN
        CALL zmeq(mb,mc)
        IF (ntrace/=0) THEN
          CALL zmntr(1,mb,mb,1)
          IF (abs(ntrace)>=1 .AND. ncall<=lvltrc) THEN
            IF (ntrace<0) THEN
              CALL zmntrj(mc,ndig)
            ELSE
              CALL zmprnt(mc)
            END IF
          END IF
        END IF
        ncall = ncall - 1
        RETURN
      END IF
      marz = ma(0)
      maiz = ma(kptimu)
      kaccsw = 0
      krsave = krad
      krad = 1

      CALL zmeq2(ma,ma,ndsave,ndig,1)

!             Check for special cases.

      IF (ma(2)==0 .AND. ma(kptimu+2)==0) THEN
        CALL zmi2m(1,mz01)
        CALL zmi2m(0,mc)
        GO TO 10
      ELSE IF (ma(kptimu+2)==0) THEN
        CALL fmchsh(ma,mz01,mc)
        CALL fmi2m(0,mz01(kptimu))
        CALL fmi2m(0,mc(kptimu))
        GO TO 10
      ELSE IF (ma(2)==0) THEN
        CALL fmcssn(ma(kptimu),mz01,mc(kptimu))
        CALL fmi2m(0,mz01(kptimu))
        CALL fmi2m(0,mc)
        GO TO 10
      END IF

!             Find SINH(REAL(MA)) and COSH(REAL(MA)).

      CALL fmchsh(ma,mz02,mz02(kptimu))

!             Find SIN(IMAG(MA)) and COS(IMAG(MA)).

      CALL fmcssn(ma(kptimu),mz03,mz03(kptimu))

!             COSH(MA) =  COSH(REAL(MA))*COS(IMAG(MA)) +
!                         SINH(REAL(MA))*SIN(IMAG(MA)) i

      CALL fmmpy(mz02,mz03,mz01)
      CALL fmmpy(mz02(kptimu),mz03(kptimu),mz01(kptimu))

!             SINH(MA) =  SINH(REAL(MA))*COS(IMAG(MA)) +
!                         COSH(REAL(MA))*SIN(IMAG(MA)) i

      CALL fmmpy(mz02(kptimu),mz03,mc)
      CALL fmmpy(mz02,mz03(kptimu),mc(kptimu))

10    maccmb = mz01(0)
      ma(0) = marz
      mz01(0) = min(maccmb,marz,maiz)
      maccmb = mz01(kptimu)
      ma(kptimu) = maiz
      mz01(kptimu) = min(maccmb,marz,maiz)
      mc(0) = mz01(0)
      mc(kptimu) = mz01(kptimu)
      kaccsw = kasave
      CALL zmeq2(mc,mc,ndig,ndsave,1)
      CALL zmexit(mz01,mb,ndsave,mxsave,kasave,kovun,0)
      IF (ntrace/=0) THEN
        IF (abs(ntrace)>=1 .AND. ncall+1<=lvltrc) THEN
          IF (ntrace<0) THEN
            CALL zmntrj(mc,ndig)
          ELSE
            CALL zmprnt(mc)
          END IF
        END IF
      END IF
      krad = krsave
      RETURN
    END SUBROUTINE zmchsh
    SUBROUTINE zmcmpx(mafm,mbfm,mc)

!  MC = COMPLEX( MAFM , MBFM )

!  MAFM and MBFM are real FM numbers, MC is a complex ZM number.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: mafm(0:lunpck), mbfm(0:lunpck), mc(0:lunpkz)
! ..
! .. External Subroutines ..
      EXTERNAL fmeq, fmntr, zmntr
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
      kflag = 0
      ncall = ncall + 1
      namest(ncall) = 'ZMCMPX'
      IF (ntrace/=0) CALL fmntr(2,mafm,mbfm,2)

      CALL fmeq(mafm,mc)
      CALL fmeq(mbfm,mc(kptimu))

      IF (ntrace/=0) CALL zmntr(1,mc,mc,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmcmpx
    SUBROUTINE zmconj(ma,mb)

!  MB = CONJG(MA)

!  Complex conjugate.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. External Subroutines ..
      EXTERNAL fmeq, zmntr
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
      kflag = 0
      ncall = ncall + 1
      namest(ncall) = 'ZMCONJ'
      IF (ntrace/=0) CALL zmntr(2,ma,ma,1)

      CALL fmeq(ma,mb)
      CALL fmeq(ma(kptimu),mb(kptimu))
      IF (mb(kptimu+1)/=munkno) mb(kptimu+2) = -mb(kptimu+2)

      IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmconj
    SUBROUTINE zmcos(ma,mb)

!  MB = COS(MA).

      IMPLICIT NONE

!             Scratch array usage during ZMCOS:   M01 - M06, MZ01 - MZ03

! .. Intrinsic Functions ..
      INTRINSIC min
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccmb, maiz, marz, mxsave
      INTEGER :: kasave, kovun, kreslt, krsave, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmchsh, fmcos, fmcosh, fmcssn, fmi2m, fmmpy, zmentr, zmeq2, &
        zmexit, zmi2m
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
        mpisav(0:lunpck), mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      CALL zmentr('ZMCOS ',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
      IF (kreslt/=0) RETURN
      marz = ma(0)
      maiz = ma(kptimu)
      kaccsw = 0
      krsave = krad
      krad = 1

      CALL zmeq2(ma,ma,ndsave,ndig,1)

!             Check for special cases.

      IF (ma(2)==0 .AND. ma(kptimu+2)==0) THEN
        CALL zmi2m(1,mz01)
        GO TO 10
      ELSE IF (ma(kptimu+2)==0) THEN
        CALL fmcos(ma,mz01)
        CALL fmi2m(0,mz01(kptimu))
        GO TO 10
      ELSE IF (ma(2)==0) THEN
        CALL fmcosh(ma(kptimu),mz01)
        CALL fmi2m(0,mz01(kptimu))
        GO TO 10
      END IF

!             Find COS(REAL(MA)) and SIN(REAL(MA)).

      CALL fmcssn(ma,mz01,mz01(kptimu))

!             Find COSH(IMAG(MA)) and SINH(IMAG(MA)).

      CALL fmchsh(ma(kptimu),m05,m06)

!             COS(MA) =  COS(REAL(MA))*COSH(IMAG(MA)) -
!                        SIN(REAL(MA))*SINH(IMAG(MA)) i

      CALL fmmpy(mz01,m05,mz01)
      IF (m06(1)/=munkno) m06(2) = -m06(2)
      CALL fmmpy(mz01(kptimu),m06,mz01(kptimu))

10    maccmb = mz01(0)
      ma(0) = marz
      mz01(0) = min(maccmb,marz,maiz)
      maccmb = mz01(kptimu)
      ma(kptimu) = maiz
      mz01(kptimu) = min(maccmb,marz,maiz)
      CALL zmexit(mz01,mb,ndsave,mxsave,kasave,kovun,0)
      krad = krsave
      RETURN
    END SUBROUTINE zmcos
    SUBROUTINE zmcosh(ma,mb)

!  MB = COSH(MA).

      IMPLICIT NONE

!             Scratch array usage during ZMCOSH:  M01 - M06, MZ01 - MZ03

! .. Intrinsic Functions ..
      INTRINSIC min
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccmb, maiz, marz, mxsave
      INTEGER :: kasave, kovun, kreslt, krsave, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmchsh, fmcos, fmcosh, fmcssn, fmi2m, fmmpy, zmentr, zmeq2, &
        zmexit, zmi2m
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
        mpisav(0:lunpck), mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      CALL zmentr('ZMCOSH',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
      IF (kreslt/=0) RETURN
      marz = ma(0)
      maiz = ma(kptimu)
      kaccsw = 0
      krsave = krad
      krad = 1

      CALL zmeq2(ma,ma,ndsave,ndig,1)

!             Check for special cases.

      IF (ma(2)==0 .AND. ma(kptimu+2)==0) THEN
        CALL zmi2m(1,mz01)
        GO TO 10
      ELSE IF (ma(2)==0) THEN
        CALL fmcos(ma(kptimu),mz01)
        CALL fmi2m(0,mz01(kptimu))
        GO TO 10
      ELSE IF (ma(kptimu+2)==0) THEN
        CALL fmcosh(ma,mz01)
        CALL fmi2m(0,mz01(kptimu))
        GO TO 10
      END IF

!             Find COS(IMAG(MA)) and SIN(IMAG(MA)).

      CALL fmcssn(ma(kptimu),mz01,mz01(kptimu))

!             Find COSH(REAL(MA)) and SINH(REAL(MA)).

      CALL fmchsh(ma,m05,m06)

!             COSH(MA) =  COSH(REAL(MA))*COS(IMAG(MA)) +
!                         SINH(REAL(MA))*SIN(IMAG(MA)) i

      CALL fmmpy(mz01,m05,mz01)
      CALL fmmpy(mz01(kptimu),m06,mz01(kptimu))

10    maccmb = mz01(0)
      ma(0) = marz
      mz01(0) = min(maccmb,marz,maiz)
      maccmb = mz01(kptimu)
      ma(kptimu) = maiz
      mz01(kptimu) = min(maccmb,marz,maiz)
      CALL zmexit(mz01,mb,ndsave,mxsave,kasave,kovun,0)
      krad = krsave
      RETURN
    END SUBROUTINE zmcosh
    SUBROUTINE zmcssn(ma,mb,mc)

!  MB = COS(MA),    MC = SIN(MA).

!  If both the sine and cosine are needed, this routine is faster
!  than calling both ZMCOS and ZMSIN.

!  MB and MC must be distinct arrays.

      IMPLICIT NONE

!             Scratch array usage during ZMCSSN:  M01 - M06, MZ01 - MZ04

! .. Intrinsic Functions ..
      INTRINSIC abs, min
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz), mc(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccmb, maiz, marz, mxsave
      INTEGER :: kasave, kovun, kreslt, krsave, ncsave, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmchsh, fmcssn, fmi2m, fmmpy, zmentr, zmeq, zmeq2, zmexit, &
        zmi2m, zmntr, zmntrj, zmprnt
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
        mpisav(0:lunpck), mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      ncsave = ncall
      CALL zmentr('ZMCSSN',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
      ncall = ncsave + 1
      IF (kreslt/=0) THEN
        CALL zmeq(mb,mc)
        IF (ntrace/=0) THEN
          CALL zmntr(1,mb,mb,1)
          IF (abs(ntrace)>=1 .AND. ncall<=lvltrc) THEN
            IF (ntrace<0) THEN
              CALL zmntrj(mc,ndig)
            ELSE
              CALL zmprnt(mc)
            END IF
          END IF
        END IF
        ncall = ncall - 1
        RETURN
      END IF
      marz = ma(0)
      maiz = ma(kptimu)
      kaccsw = 0
      krsave = krad
      krad = 1

      CALL zmeq2(ma,ma,ndsave,ndig,1)

!             Check for special cases.

      IF (ma(2)==0 .AND. ma(kptimu+2)==0) THEN
        CALL zmi2m(1,mz01)
        CALL zmi2m(0,mc)
        GO TO 10
      ELSE IF (ma(kptimu+2)==0) THEN
        CALL fmcssn(ma,mz01,mc)
        CALL fmi2m(0,mz01(kptimu))
        CALL fmi2m(0,mc(kptimu))
        GO TO 10
      ELSE IF (ma(2)==0) THEN
        CALL fmchsh(ma(kptimu),mz01,mc(kptimu))
        CALL fmi2m(0,mz01(kptimu))
        CALL fmi2m(0,mc)
        GO TO 10
      END IF

!             Find SIN(REAL(MA)) and COS(REAL(MA)).

      CALL fmcssn(ma,mz02,mz02(kptimu))

!             Find SINH(IMAG(MA)) and COSH(IMAG(MA)).

      CALL fmchsh(ma(kptimu),mz03,mz03(kptimu))

!             COS(MA) =  COS(REAL(MA))*COSH(IMAG(MA)) -
!                        SIN(REAL(MA))*SINH(IMAG(MA)) i

      CALL fmmpy(mz02,mz03,mz01)
      CALL fmmpy(mz02(kptimu),mz03(kptimu),mz01(kptimu))
      IF (mz01(kptimu+1)/=munkno) mz01(kptimu+2) = -mz01(kptimu+2)

!             SIN(MA) =  SIN(REAL(MA))*COSH(IMAG(MA)) +
!                        COS(REAL(MA))*SINH(IMAG(MA)) i

      CALL fmmpy(mz02(kptimu),mz03,mc)
      CALL fmmpy(mz02,mz03(kptimu),mc(kptimu))

10    maccmb = mz01(0)
      ma(0) = marz
      mz01(0) = min(maccmb,marz,maiz)
      maccmb = mz01(kptimu)
      ma(kptimu) = maiz
      mz01(kptimu) = min(maccmb,marz,maiz)
      mc(0) = mz01(0)
      mc(kptimu) = mz01(kptimu)
      kaccsw = kasave
      CALL zmeq2(mc,mc,ndig,ndsave,1)
      CALL zmexit(mz01,mb,ndsave,mxsave,kasave,kovun,0)
      IF (ntrace/=0) THEN
        IF (abs(ntrace)>=1 .AND. ncall+1<=lvltrc) THEN
          IF (ntrace<0) THEN
            CALL zmntrj(mc,ndig)
          ELSE
            CALL zmprnt(mc)
          END IF
        END IF
      END IF
      krad = krsave
      RETURN
    END SUBROUTINE zmcssn
    SUBROUTINE zmdiv(ma,mb,mc)

!  MC = MA / MB

      IMPLICIT NONE

!             Scratch array usage during ZMDIV:   M01 - M04, MZ01

! .. Intrinsic Functions ..
      INTRINSIC abs, int, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz), mc(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccmb, maiz, marz, mbiz, mbrz, mxsave, mz11sv, mz1ksv, &
        mzero
      INTEGER :: iextra, j, kasave, kovun, kreslt, kwrnsv, ndgsv2, ndsave, &
        ngoal, ntrsav
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmdiv, fmdivd, fmmpye, fmsub, zmentr, zmeq2, zmi2m, &
        zmntr, zmrslt, zmwarn
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
        mpisav(0:lunpck), mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      IF (abs(ma(1))>mexpab .OR. abs(ma(kptimu+1))>mexpab .OR. abs(mb(1))> &
          mexpab .OR. abs(mb(kptimu+1))>mexpab .OR. kdebug>=1) THEN
        CALL zmentr('ZMDIV ',ma,mb,2,mc,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        IF (ntrace/=0) THEN
          namest(ncall) = 'ZMDIV '
          CALL zmntr(2,ma,mb,2)
        END IF
        ndsave = ndig
        IF (ncall==1) THEN
          ndig = max(ndig+ngrd52,2)
          IF (ndig>ndg2mx) THEN
            namest(ncall) = 'ZMDIV '
            kflag = -9
            CALL zmwarn
            kreslt = 12
            ndig = ndsave
            CALL zmrslt(mc,kreslt)
            IF (ntrace/=0) CALL zmntr(1,mc,mc,1)
            ncall = ncall - 1
            RETURN
          END IF
          IF (mbase>=100*abs(ma(2)) .OR. mbase>=100*abs(ma(kptimu+2))) THEN
            ndig = min(ndig+1,ndg2mx)
          ELSE IF (mbase>=100*abs(mb(2)) .OR. mbase>=100*abs(mb(kptimu+ &
              2))) THEN
            ndig = min(ndig+1,ndg2mx)
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 1
        mxsave = mxexp
        mxexp = mxexp2
        kovun = 0
      END IF

      marz = ma(0)
      mbrz = mb(0)
      maiz = ma(kptimu)
      mbiz = mb(kptimu)
      mzero = 0
      ntrsav = ntrace
      ntrace = 0
      kwrnsv = kwarn
      kwarn = 0
      iextra = 0
      mz11sv = -munkno
      mz1ksv = -munkno

10    DO 20 j = ndsave + 2, ndig + 1
        ma(j) = mzero
        mb(j) = mzero
        ma(kptimu+j) = mzero
        mb(kptimu+j) = mzero
20    CONTINUE
      IF (ncall==1) THEN
        ma(0) = nint(ndig*alogm2)
        mb(0) = ma(0)
        ma(kptimu) = ma(0)
        mb(kptimu) = ma(0)
      END IF

!             Check for special cases.

      IF (mb(kptimu+2)==0) THEN
        CALL fmdivd(ma,ma(kptimu),mb,mz01,mz01(kptimu))
        GO TO 70
      ELSE IF (mb(2)==0) THEN
        CALL fmdivd(ma(kptimu),ma,mb(kptimu),mz01,mz01(kptimu))
        IF (mz01(kptimu+1)/=munkno) mz01(kptimu+2) = -mz01(kptimu+2)
        GO TO 70
      END IF
      IF (ma(1)==mb(1) .AND. ma(2)==mb(2)) THEN
        IF (ma(kptimu+1)==mb(kptimu+1) .AND. ma(kptimu+2)==mb(kptimu+2)) THEN
          DO 30 j = 3, ndsave + 1
            IF (ma(j)/=mb(j)) GO TO 50
            IF (ma(kptimu+j)/=mb(kptimu+j)) GO TO 50
30        CONTINUE
          IF (abs(ma(1))<mexpov .AND. abs(ma(kptimu+ &
              1))<mexpov .AND. abs(mb(1))<mexpov .AND. abs(mb(kptimu+ &
              1))<mexpov) THEN
            CALL zmi2m(1,mz01)
            GO TO 70
          END IF
        END IF
      END IF
      IF (ma(1)==mb(1) .AND. (-ma(2))==mb(2)) THEN
        IF (ma(kptimu+1)==mb(kptimu+1) .AND. (-ma(kptimu+ &
            2))==mb(kptimu+2)) THEN
          DO 40 j = 3, ndsave + 1
            IF (ma(j)/=mb(j)) GO TO 50
            IF (ma(kptimu+j)/=mb(kptimu+j)) GO TO 50
40        CONTINUE
          IF (abs(ma(1))<mexpov .AND. abs(ma(kptimu+ &
              1))<mexpov .AND. abs(mb(1))<mexpov .AND. abs(mb(kptimu+ &
              1))<mexpov) THEN
            CALL zmi2m(-1,mz01)
            GO TO 70
          END IF
        END IF
      END IF
50    IF (mz11sv/=-munkno) THEN

!             If a retry is being done due to cancellation, try a slower
!             but more stable form of the division formula.

        CALL fmmpye(mb,ma,ma(kptimu),mb,mz01,mz01(kptimu),m03)
        CALL fmmpye(mb(kptimu),ma(kptimu),ma,mb(kptimu),m01,m02,m04)
        CALL fmadd(m03,m04,m04)
        CALL fmadd(mz01,m01,mz01)
        CALL fmsub(mz01(kptimu),m02,mz01(kptimu))
        CALL fmdivd(mz01,mz01(kptimu),m04,mz01,mz01(kptimu))
        IF (abs(mz01(1))<mexpov .AND. abs(mz01(kptimu+1))<mexpov) GO TO 70
      END IF

!             Normal method for  ( a + b i ) / ( c + d i ):

!             If  abs(c) << abs(d)  Then

!                 P = c / d
!                 result = ( a*P + b )/( c*P + d ) +
!                          ( b*P - a )/( c*P + d ) i

!             Else

!                 P = d / c
!                 result = ( b*P + a )/( d*P + c ) +
!                          ( b - a*P )/( d*P + c ) i

      kaccsw = 0
      IF (mb(1)<=mb(kptimu+1)) THEN
        CALL fmdiv(mb,mb(kptimu),m04)
        CALL fmmpye(m04,ma,ma(kptimu),mb,mz01,mz01(kptimu),m03)
        IF (ma(kptimu+2)*mz01(2)<0) THEN
          kaccsw = 1
        ELSE
          kaccsw = 0
        END IF
        CALL fmadd(ma(kptimu),mz01,mz01)
        IF (m03(2)*mb(kptimu+2)<0) THEN
          kaccsw = 1
        ELSE
          kaccsw = 0
        END IF
        CALL fmadd(m03,mb(kptimu),m03)
        IF (mz01(kptimu+2)*ma(2)<0) THEN
          kaccsw = 0
        ELSE
          kaccsw = 1
        END IF
        CALL fmsub(mz01(kptimu),ma,mz01(kptimu))
        kaccsw = 0
        CALL fmdivd(mz01,mz01(kptimu),m03,mz01,mz01(kptimu))
      ELSE
        CALL fmdiv(mb(kptimu),mb,m04)
        CALL fmmpye(m04,ma(kptimu),ma,mb(kptimu),mz01,mz01(kptimu),m03)
        IF (ma(2)*mz01(2)<0) THEN
          kaccsw = 1
        ELSE
          kaccsw = 0
        END IF
        CALL fmadd(ma,mz01,mz01)
        IF (m03(2)*mb(2)<0) THEN
          kaccsw = 1
        ELSE
          kaccsw = 0
        END IF
        CALL fmadd(m03,mb,m03)
        IF (mz01(kptimu+2)*ma(kptimu+2)<0) THEN
          kaccsw = 0
        ELSE
          kaccsw = 1
        END IF
        CALL fmsub(ma(kptimu),mz01(kptimu),mz01(kptimu))
        kaccsw = 0
        CALL fmdivd(mz01,mz01(kptimu),m03,mz01,mz01(kptimu))
      END IF
      kaccsw = 1

!             Check for too much cancellation.

      IF (ncall<=1) THEN
        ngoal = int(real(ndsave)*alogm2) + 7
      ELSE
        ngoal = int(-mxexp2)
      END IF
      IF (mz01(0)<=ngoal .OR. mz01(kptimu)<=ngoal) THEN
        IF (mz11sv-mz01(1)>=iextra-1 .AND. mz01(kptimu)>ngoal) GO TO 70
        IF (mz1ksv-mz01(kptimu+1)>=iextra-1 .AND. mz01(0)>ngoal) GO TO 70
        IF (mz11sv>-munkno .AND. mz01(0)>ngoal .AND. mz01(kptimu+2)==0) &
          GO TO 70
        IF (mz11sv>-munkno .AND. mz01(kptimu)>ngoal .AND. mz01(2)==0) GO TO 70
        iextra = int(real(max(ngoal-mz01(0),ngoal-mz01(kptimu)))/alogm2+23.03/ &
          alogmb) + 1
        mz11sv = mz01(1)
        mz1ksv = mz01(kptimu+1)
        ndig = ndig + iextra
        IF (ndig>ndg2mx) THEN
          namest(ncall) = 'ZMDIV '
          kflag = -9
          CALL zmwarn
          mz01(1) = munkno
          mz01(2) = 1
          mz01(kptimu+1) = munkno
          mz01(kptimu+2) = 1
          DO 60 j = 2, ndsave
            mz01(j+1) = 0
            mz01(kptimu+j+1) = 0
60        CONTINUE
          ndig = ndig - iextra
          mz01(0) = nint(ndig*alogm2)
          mz01(kptimu) = nint(ndig*alogm2)
          GO TO 70
        END IF
        GO TO 10
      END IF

70    mxexp = mxsave
      ntrace = ntrsav
      ndgsv2 = ndig
      ndig = ndsave
      kwarn = kwrnsv
      maccmb = mz01(0)
      ma(0) = marz
      mb(0) = mbrz
      mz01(0) = min(maccmb,marz,maiz,mbrz,mbiz)
      maccmb = mz01(kptimu)
      ma(kptimu) = maiz
      mb(kptimu) = mbiz
      mz01(kptimu) = min(maccmb,marz,maiz,mbrz,mbiz)
      CALL zmeq2(mz01,mc,ndgsv2,ndsave,0)
      IF (mc(1)>=mexpov .OR. mc(1)<=-mexpov .OR. mc(kptimu+1)>=mexpov .OR. &
          mc(kptimu+1)<=-mexpov) THEN
        IF (mc(1)==munkno .OR. mc(kptimu+1)==munkno) THEN
          kflag = -4
        ELSE IF (mc(1)==mexpov .OR. mc(kptimu+1)==mexpov) THEN
          kflag = -5
        ELSE IF (mc(1)==mexpun .OR. mc(kptimu+1)==mexpun) THEN
          kflag = -6
        END IF
        IF ((mc(1)==munkno) .OR. (mc(kptimu+1)==munkno) .OR. (mc(1)==mexpun &
            .AND. kovun==0) .OR. (mc(kptimu+1)==mexpun .AND. kovun==0) .OR. ( &
            mc(1)==mexpov .AND. kovun==0) .OR. (mc(kptimu+ &
            1)==mexpov .AND. kovun==0)) THEN
          namest(ncall) = 'ZMDIV '
          CALL zmwarn
        END IF
      END IF
      IF (ntrace/=0) CALL zmntr(1,mc,mc,1)
      kaccsw = kasave
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmdiv
    SUBROUTINE zmdivi(ma,integ,mb)

!  MB = MA / INTEG        Divide by one-word (real) integer.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: integ
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: mxsave
      INTEGER :: kasave, kovun, kreslt, kwrnsv, ndsave, ntrsav
! ..
! .. External Subroutines ..
      EXTERNAL fmdivi, fmntri, zmentr, zmntr, zmwarn
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
      IF (abs(ma(1))>mexpab .OR. abs(ma(kptimu+1))>mexpab .OR. kdebug>=1) THEN
        CALL zmentr('ZMDIVI',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
        ndig = ndsave
        mxexp = mxsave
        kaccsw = kasave
        ntrsav = ntrace
      ELSE
        ncall = ncall + 1
        IF (ntrace/=0) THEN
          namest(ncall) = 'ZMDIVI'
          CALL zmntr(2,ma,ma,1)
          CALL fmntri(2,integ,0)
        END IF
        kovun = 0
      END IF

!             Force FMDIVI to use more guard digits for user calls.

      ncall = ncall - 1
      ntrsav = ntrace
      ntrace = 0
      kwrnsv = kwarn
      kwarn = 0

      CALL fmdivi(ma,integ,mb)
      CALL fmdivi(ma(kptimu),integ,mb(kptimu))

      ntrace = ntrsav
      kwarn = kwrnsv
      ncall = ncall + 1
      IF (ntrace/=0) namest(ncall) = 'ZMDIVI'
      IF (mb(1)==munkno .OR. mb(kptimu+1)==munkno) THEN
        kflag = -4
      ELSE IF (mb(1)==mexpov .OR. mb(kptimu+1)==mexpov) THEN
        kflag = -5
      ELSE IF (mb(1)==mexpun .OR. mb(kptimu+1)==mexpun) THEN
        kflag = -6
      END IF
      IF ((mb(1)==munkno) .OR. (mb(kptimu+1)==munkno) .OR. (mb(1)==mexpun &
          .AND. kovun==0) .OR. (mb(kptimu+1)==mexpun .AND. kovun==0) .OR. (mb( &
          1)==mexpov .AND. kovun==0) .OR. (mb(kptimu+ &
          1)==mexpov .AND. kovun==0)) THEN
        namest(ncall) = 'ZMDIVI'
        CALL zmwarn
      END IF
      IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmdivi
    SUBROUTINE zmentr(nroutn,ma,mb,nargs,mc,kreslt,ndsave,mxsave,kasave,kovun)

!  Do the argument checking and increasing of precision, overflow
!  threshold, etc., upon entry to a ZM routine.

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
      INTRINSIC abs, int, max, min, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL (KIND(0.0D0)) :: mxsave
      INTEGER :: kasave, kovun, kreslt, nargs, ndsave
      CHARACTER (6) :: nroutn
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz), mc(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: mbs
      INTEGER :: j, kwrnsv, nds
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, zmi2m, zmntr, zmrslt, zmwarn
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
      kreslt = 0
      ncall = ncall + 1
      kflag = 0
      namest(ncall) = nroutn
      IF (ntrace/=0) CALL zmntr(2,ma,mb,nargs)

      IF (mblogs/=mbase) CALL fmcons
      kovun = 0
      IF (ma(1)==mexpov .OR. ma(1)==mexpun .OR. ma(kptimu+1)==mexpov .OR. &
        ma(kptimu+1)==mexpun) kovun = 1
      IF (nargs==2) THEN
        IF (mb(1)==mexpov .OR. mb(1)==mexpun .OR. mb(kptimu+1)==mexpov .OR. &
          mb(kptimu+1)==mexpun) kovun = 1
      END IF
      kasave = kaccsw
      mxsave = mxexp

!             Check the validity of parameters if this is a user call.

      IF (ncall>1 .AND. kdebug==0) GO TO 70

!             Check NDIG.

      IF (ndig<2 .OR. ndig>ndigmx) THEN
        kflag = -1
        CALL zmwarn
        nds = ndig
        IF (ndig<2) ndig = 2
        IF (ndig>ndigmx) ndig = ndigmx
        WRITE (kw,90000) nds, ndig
        kreslt = 12
        GO TO 70
      END IF

!             Check MBASE.

      IF (mbase<2 .OR. mbase>mxbase) THEN
        kflag = -2
        CALL zmwarn
        mbs = mbase
        IF (mbase<2) mbase = 2
        IF (mbase>mxbase) mbase = mxbase
        WRITE (kw,90010) int(mbs), int(mbase)
        CALL fmcons
        kreslt = 12
        GO TO 70
      END IF

!             Check exponent range.

      IF (ma(1)>mxexp+1 .OR. ma(1)<-mxexp) THEN
        IF ((abs(ma(1))/=mexpov .AND. abs(ma(1))/=munkno) .OR. abs(ma(2))/=1) &
            THEN
          kflag = -3
          CALL zmwarn
          CALL zmi2m(0,ma)
          ma(1) = munkno
          ma(2) = 1
          ma(kptimu+1) = munkno
          ma(kptimu+2) = 1
          ma(0) = nint(ndig*alogm2)
          ma(kptimu) = nint(ndig*alogm2)
          kreslt = 12
          GO TO 70
        END IF
      END IF
      IF (ma(kptimu+1)>mxexp+1 .OR. ma(kptimu+1)<-mxexp) THEN
        IF ((abs(ma(kptimu+1))/=mexpov .AND. abs(ma(kptimu+1))/=munkno) .OR. &
            abs(ma(kptimu+2))/=1) THEN
          kflag = -3
          CALL zmwarn
          CALL zmi2m(0,ma)
          ma(1) = munkno
          ma(2) = 1
          ma(kptimu+1) = munkno
          ma(kptimu+2) = 1
          ma(0) = nint(ndig*alogm2)
          ma(kptimu) = nint(ndig*alogm2)
          kreslt = 12
          GO TO 70
        END IF
      END IF
      IF (nargs==2) THEN
        IF (mb(1)>mxexp+1 .OR. mb(1)<-mxexp) THEN
          IF ((abs(mb(1))/=mexpov .AND. abs(mb(1))/=munkno) .OR. abs(mb( &
              2))/=1) THEN
            kflag = -3
            CALL zmwarn
            CALL zmi2m(0,mb)
            mb(1) = munkno
            mb(2) = 1
            mb(kptimu+1) = munkno
            mb(kptimu+2) = 1
            mb(0) = nint(ndig*alogm2)
            mb(kptimu) = nint(ndig*alogm2)
            kreslt = 12
            GO TO 70
          END IF
        END IF
        IF (mb(kptimu+1)>mxexp+1 .OR. mb(kptimu+1)<-mxexp) THEN
          IF ((abs(mb(kptimu+1))/=mexpov .AND. abs(mb(kptimu+1))/=munkno) .OR. &
              abs(mb(kptimu+2))/=1) THEN
            kflag = -3
            CALL zmwarn
            CALL zmi2m(0,mb)
            mb(1) = munkno
            mb(2) = 1
            mb(kptimu+1) = munkno
            mb(kptimu+2) = 1
            mb(0) = nint(ndig*alogm2)
            mb(kptimu) = nint(ndig*alogm2)
            kreslt = 12
            GO TO 70
          END IF
        END IF
      END IF

!             Check for properly normalized digits in the
!             input arguments.

      IF (abs(ma(1)-int(ma(1)))/=0) kflag = 1
      IF (abs(ma(kptimu+1)-int(ma(kptimu+1)))/=0) kflag = kptimu + 1
      IF (ma(2)<=(-mbase) .OR. ma(2)>=mbase .OR. abs(ma(2)-int(ma(2)))/=0) &
        kflag = 2
      IF (ma(kptimu+2)<=(-mbase) .OR. ma(kptimu+2)>=mbase .OR. abs(ma(kptimu+ &
        2)-int(ma(kptimu+2)))/=0) kflag = kptimu + 2
      IF (kdebug==0) GO TO 30
      DO 10 j = 3, ndig + 1
        IF (ma(j)<0 .OR. ma(j)>=mbase .OR. abs(ma(j)-int(ma(j)))/=0) THEN
          kflag = j
          GO TO 30
        END IF
10    CONTINUE
      DO 20 j = kptimu + 3, kptimu + ndig + 1
        IF (ma(j)<0 .OR. ma(j)>=mbase .OR. abs(ma(j)-int(ma(j)))/=0) THEN
          kflag = j
          GO TO 30
        END IF
20    CONTINUE
30    IF (kflag/=0) THEN
        j = kflag
        mbs = ma(j)
        CALL zmi2m(0,ma)
        kflag = -4
        kwrnsv = kwarn
        IF (kwarn>=2) kwarn = 1
        CALL zmwarn
        kwarn = kwrnsv
        IF (kwarn>=1) THEN
          IF (j<kptimu) THEN
            WRITE (kw,*) ' First invalid array element:  MA(', j, ') = ', mbs
          ELSE
            WRITE (kw,*) ' First invalid array element:  MA(', kptimu, '+', &
              j - kptimu, ') = ', mbs
          END IF
        END IF
        ma(1) = munkno
        ma(2) = 1
        ma(kptimu+1) = munkno
        ma(kptimu+2) = 1
        ma(0) = nint(ndig*alogm2)
        ma(kptimu) = nint(ndig*alogm2)
        IF (kwarn>=2) THEN
          STOP
        END IF
        kreslt = 12
        GO TO 70
      END IF
      IF (nargs==2) THEN
        IF (abs(mb(1)-int(mb(1)))/=0) kflag = 1
        IF (abs(mb(kptimu+1)-int(mb(kptimu+1)))/=0) kflag = kptimu + 1
        IF (mb(2)<=(-mbase) .OR. mb(2)>=mbase .OR. abs(mb(2)-int(mb(2)))/=0) &
          kflag = 2
        IF (mb(kptimu+2)<=(-mbase) .OR. mb(kptimu+2)>=mbase .OR. abs(mb( &
          kptimu+2)-int(mb(kptimu+2)))/=0) kflag = kptimu + 2
        IF (kdebug==0) GO TO 60
        DO 40 j = 3, ndig + 1
          IF (mb(j)<0 .OR. mb(j)>=mbase .OR. abs(mb(j)-int(mb(j)))/=0) THEN
            kflag = j
            GO TO 60
          END IF
40      CONTINUE
        DO 50 j = kptimu + 3, kptimu + ndig + 1
          IF (mb(j)<0 .OR. mb(j)>=mbase .OR. abs(mb(j)-int(mb(j)))/=0) THEN
            kflag = j
            GO TO 60
          END IF
50      CONTINUE
60      IF (kflag/=0) THEN
          j = kflag
          mbs = mb(j)
          CALL zmi2m(0,mb)
          kflag = -4
          kwrnsv = kwarn
          IF (kwarn>=2) kwarn = 1
          CALL zmwarn
          kwarn = kwrnsv
          IF (kwarn>=1) THEN
            IF (j<kptimu) THEN
              WRITE (kw,*) ' First invalid array element:  MB(', j, ') = ', &
                mbs
            ELSE
              WRITE (kw,*) ' First invalid array element:  MB(', kptimu, '+', &
                j - kptimu, ') = ', mbs
            END IF
          END IF
          mb(1) = munkno
          mb(2) = 1
          mb(kptimu+1) = munkno
          mb(kptimu+2) = 1
          mb(0) = nint(ndig*alogm2)
          mb(kptimu) = nint(ndig*alogm2)
          IF (kwarn>=2) THEN
            STOP
          END IF
          kreslt = 12
          GO TO 70
        END IF
      END IF

!             Increase the working precision.

70    ndsave = ndig
      IF (ncall==1) THEN
        ndig = max(ndig+ngrd52,2)
        IF (ndig>ndg2mx) THEN
          kflag = -9
          CALL zmwarn
          kreslt = 12
          ndig = ndsave
        END IF
        IF (mbase>=100*abs(ma(2)) .OR. mbase>=100*abs(ma(kptimu+2))) THEN
          ndig = min(ndig+1,ndg2mx)
        ELSE IF (nargs==2 .AND. (mbase>=100*abs(mb(2)) .OR. mbase>=100*abs( &
            mb(kptimu+2)))) THEN
          ndig = min(ndig+1,ndg2mx)
        END IF
      END IF
      IF ((ma(1)==munkno .AND. ma(kptimu+1)==munkno) .OR. (mb(1)==munkno .AND. &
          mb(kptimu+1)==munkno)) THEN
        kflag = -4
        kreslt = 12
      END IF
      IF (kreslt/=0) THEN
        ndig = ndsave
        CALL zmrslt(mc,kreslt)
        IF (ntrace/=0) CALL zmntr(1,mc,mc,1)
        ncall = ncall - 1
        RETURN
      END IF

      kaccsw = 1

!             Extend the overflow/underflow threshold.

      mxexp = mxexp2
      RETURN
90000 FORMAT (' NDIG was',I10,'.  It has been changed to',I10,'.')
90010 FORMAT (' MBASE was',I10,'.  It has been changed to',I10,'.')
    END SUBROUTINE zmentr
    SUBROUTINE zmeq(ma,mb)

!  MB = MA

!  This is the standard form of equality, where MA and MB both
!  have precision NDIG.  Use ZMEQU for assignments that also
!  change precision.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. External Subroutines ..
      EXTERNAL fmeq
! ..
      CALL fmeq(ma,mb)
      CALL fmeq(ma(kptimu),mb(kptimu))
      RETURN
    END SUBROUTINE zmeq
    SUBROUTINE zmeq2(ma,mb,nda,ndb,ksame)

!  Set MB (having NDB digits) equal to MA (having NDA digits).

!  If MA and MB are the same array, setting KSAME = 1 before calling
!  ZMEQ2 gives faster performance.

!  If MB has less precision than MA, the result is rounded to
!  NDB digits.

!  If MB has more precision, the result has zero digits padded on the
!  right.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ksame, nda, ndb
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. External Subroutines ..
      EXTERNAL fmeq2
! ..
      CALL fmeq2(ma,mb,nda,ndb,ksame)
      CALL fmeq2(ma(kptimu),mb(kptimu),nda,ndb,ksame)
      RETURN
    END SUBROUTINE zmeq2
    SUBROUTINE zmequ(ma,mb,nda,ndb)

!  Set MB (having NDB digits) equal to MA (having NDA digits).

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: nda, ndb
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. External Subroutines ..
      EXTERNAL fmeq2
! ..
      CALL fmeq2(ma,mb,nda,ndb,0)
      CALL fmeq2(ma(kptimu),mb(kptimu),nda,ndb,0)
      RETURN
    END SUBROUTINE zmequ
    SUBROUTINE zmexit(mt,mc,ndsave,mxsave,kasave,kovun,ksame)

!  Upon exit from an ZM routine the result MT (having precision NDIG)
!  is rounded and returned in MC (having precision NDSAVE).
!  The values of NDIG, MXEXP, and KACCSW are restored to the values
!  NDSAVE,MXSAVE,KASAVE.
!  KSAME is 1 if MT and MC are the same array in the calling routine.
!  KOVUN is nonzero if one of the routine's input arguments was overflow
!  or underflow.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL (KIND(0.0D0)) :: mxsave
      INTEGER :: kasave, kovun, ksame, ndsave
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: mc(0:lunpkz), mt(0:lunpkz)
! ..
! .. Local Scalars ..
      INTEGER :: kfsave, kwrnsv
! ..
! .. External Subroutines ..
      EXTERNAL zmeq2, zmntr, zmwarn
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
      kwrnsv = kwarn
      kwarn = 0
      mxexp = mxsave
      kfsave = kflag
      kaccsw = kasave
      CALL zmeq2(mt,mc,ndig,ndsave,ksame)
      IF (kflag/=-5 .AND. kflag/=-6) kflag = kfsave
      ndig = ndsave
      kwarn = kwrnsv
      IF (kflag==1) kflag = 0
      IF (mc(1)==mexpun .OR. mc(kptimu+1)==mexpun) kflag = -6
      IF (mc(1)==mexpov .OR. mc(kptimu+1)==mexpov) kflag = -5
      IF (mc(1)==munkno .OR. mc(kptimu+1)==munkno) THEN
        IF (kflag/=-9) kflag = -4
      END IF
      IF ((mc(1)==munkno .AND. kflag/=-9) .OR. (mc(kptimu+ &
        1)==munkno .AND. kflag/=-9) .OR. (mc(1)==mexpun .AND. kovun==0) .OR. ( &
        mc(kptimu+1)==mexpun .AND. kovun==0) .OR. (mc(1)==mexpov .AND. kovun== &
        0) .OR. (mc(kptimu+1)==mexpov .AND. kovun==0)) CALL zmwarn
      IF (ntrace/=0) CALL zmntr(1,mc,mc,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmexit
    SUBROUTINE zmexi2(mtfm,mcfm,ndsave,mxsave,kasave,kovun,ksame)

!  This routine is used upon exit for complex functions that
!  return real FM results.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      REAL (KIND(0.0D0)) :: mxsave
      INTEGER :: kasave, kovun, ksame, ndsave
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: mcfm(0:lunpck), mtfm(0:lunpck)
! ..
! .. Local Scalars ..
      INTEGER :: kfsave, kwrnsv
! ..
! .. External Subroutines ..
      EXTERNAL fmeq2, zmntr2, zmwarn
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
      kwrnsv = kwarn
      kwarn = 0
      mxexp = mxsave
      kfsave = kflag
      kaccsw = kasave
      CALL fmeq2(mtfm,mcfm,ndig,ndsave,ksame)
      IF (kflag/=-5 .AND. kflag/=-6) kflag = kfsave
      ndig = ndsave
      kwarn = kwrnsv
      IF (kflag==1) kflag = 0
      IF (mcfm(1)==munkno) THEN
        IF (kflag>=0) kflag = -4
      ELSE IF (mcfm(1)==mexpov) THEN
        kflag = -5
      ELSE IF (mcfm(1)==mexpun) THEN
        kflag = -6
      END IF
      IF ((mcfm(1)==munkno .AND. kflag/=-9) .OR. (mcfm( &
        1)==mexpun .AND. kovun==0) .OR. (mcfm(1)==mexpov .AND. kovun==0)) &
        CALL zmwarn
      IF (ntrace/=0) CALL zmntr2(1,mcfm,mcfm,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmexi2
    SUBROUTINE zmexp(ma,mb)

!  MB = EXP(MA).

      IMPLICIT NONE

!             Scratch array usage during ZMEXP:   M01 - M06, MZ01

! .. Intrinsic Functions ..
      INTRINSIC min
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccmb, maiz, marz, mxsave
      INTEGER :: kasave, kovun, kreslt, krsave, kwrnsv, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmcssn, fmexp, fmi2m, fmmpyd, zmentr, zmeq2, zmexit, zmi2m
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
        mpisav(0:lunpck), mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      CALL zmentr('ZMEXP ',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
      IF (kreslt/=0) RETURN
      marz = ma(0)
      maiz = ma(kptimu)
      kaccsw = 0
      krsave = krad
      krad = 1

      CALL zmeq2(ma,ma,ndsave,ndig,1)

!             Check for special cases.

      IF (ma(2)==0 .AND. ma(kptimu+2)==0) THEN
        CALL zmi2m(1,mz01)
        GO TO 10
      ELSE IF (ma(2)==0) THEN
        CALL fmi2m(1,m06)
      ELSE
        CALL fmexp(ma,m06)
      END IF

      CALL fmcssn(ma(kptimu),mz01,mz01(kptimu))

      kwrnsv = kwarn
      kwarn = 0
      CALL fmmpyd(m06,mz01,mz01(kptimu),mz01,mz01(kptimu))
      kwarn = kwrnsv

10    maccmb = mz01(0)
      ma(0) = marz
      mz01(0) = min(maccmb,marz,maiz)
      maccmb = mz01(kptimu)
      ma(kptimu) = maiz
      mz01(kptimu) = min(maccmb,marz,maiz)
      CALL zmexit(mz01,mb,ndsave,mxsave,kasave,kovun,0)
      krad = krsave
      RETURN
    END SUBROUTINE zmexp
    SUBROUTINE zmform(form1,form2,ma,string)

!  Convert MA to STRING using FORM1 format for the real part and
!  FORM2 format for the imaginary part.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC len
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (*) :: form1, form2, string
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz)
! ..
! .. Local Scalars ..
      INTEGER :: j, kwidim, kwidre, last, lsign
! ..
! .. External Subroutines ..
      EXTERNAL fmeq, zmfpcm
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, jformz, jprntz, kaccsw, &
        kdebug, keswch, kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, &
        ncall, ndg2mx, ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff), cmbufz(lmbufz)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
      COMMON /zmbuff/cmbufz
      COMMON /zmuser/jformz, jprntz
! ..
      ncall = ncall + 1
      namest(ncall) = 'ZMFORM'
      string = ' '
      CALL zmfpcm(form1,ma,kwidre,cmbufz)
      CALL fmeq(ma(kptimu),m02)
      IF (m02(2)>=0) THEN
        lsign = 1
      ELSE
        lsign = -1
        m02(2) = -m02(2)
      END IF
      CALL zmfpcm(form2,m02,kwidim,cmbuff)

      cmbufz(kwidre+1) = ' '
      IF (lsign==1) THEN
        cmbufz(kwidre+2) = '+'
      ELSE
        cmbufz(kwidre+2) = '-'
      END IF
      cmbufz(kwidre+3) = ' '
      DO 10 j = 1, kwidim
        cmbufz(kwidre+3+j) = cmbuff(j)
10    CONTINUE
      cmbufz(kwidre+4+kwidim) = ' '
      cmbufz(kwidre+5+kwidim) = 'i'
      IF (jformz==2) cmbufz(kwidre+5+kwidim) = 'I'
      last = kwidre + kwidim + 5

      IF (last<=len(string)) THEN
        DO 20 j = 1, last
          string(j:j) = cmbufz(j)
20      CONTINUE
      ELSE
        DO 30 j = 1, last
          string(j:j) = '*'
30      CONTINUE
      END IF
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmform
    SUBROUTINE zmfpcm(form,ma,kwi,cmb)

!  Internal routine to convert MA to base 10 using FORM format.
!  The result is returned in CMB and the field width is KWI.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, dble, index, int, len, log10, max, min, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kwi
      CHARACTER (*) :: form
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(*)
      CHARACTER (1) :: cmb(lmbuff)
! ..
! .. Local Scalars ..
      INTEGER :: j, jf1sav, jf2sav, jpt, k1, k2, k3, kd, ksave, kwd, last, lb, &
        lengfm, lfirst, nd, nexp
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
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mwa(lmwa)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
! ..
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
        kwd = kwi + 11
        CALL fmnint(ma,m03)
        IF (m03(2)/=0) THEN
          CALL fmout(m03,cmb,kwd)
        ELSE
          DO 10 j = 1, kwd
            cmb(j) = ' '
10        CONTINUE
          cmb(2) = '0'
        END IF
        lfirst = 1
        last = 1
        DO 20 j = 1, kwd
          IF (cmb(kwd+1-j)/=' ') lfirst = kwd + 1 - j
          IF (cmb(j)/=' ') last = j
20      CONTINUE
        jpt = 1
        IF (last-lfirst+1>kwi) GO TO 110
        IF (last<=kwi) THEN
          DO 30 j = last, lfirst, -1
            jpt = kwi - last + j
            cmb(jpt) = cmb(j)
30        CONTINUE
          DO 40 j = 1, jpt - 1
            cmb(j) = ' '
40        CONTINUE
        ELSE
          DO 50 j = lfirst, last
            jpt = kwi - last + j
            cmb(jpt) = cmb(j)
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
        CALL fmout(ma,cmb,kwd)
        lfirst = 1
        last = 1
        DO 60 j = 1, kwd
          IF (cmb(kwd+1-j)/=' ') lfirst = kwd + 1 - j
          IF (cmb(j)/=' ') last = j
60      CONTINUE
        IF (last-lfirst+1>kwi) THEN

!             Not enough room for this F format, or FMOUT converted
!             it to E format to avoid showing no significant digits.
!             See if a shortened form will fit in E format.

          nexp = int(log10((abs(ma(1))+1)*log10(dble(mbase))+1)+1)
          nd = kwi - nexp - 5
          IF (nd<1) THEN
            GO TO 110
          ELSE
            jform1 = 0
            jform2 = nd
            CALL fmout(ma,cmb,kwi)
            lfirst = 1
            last = 1
            DO 70 j = 1, kwi
              IF (cmb(kwi+1-j)/=' ') lfirst = kwi + 1 - j
              IF (cmb(j)/=' ') last = j
70          CONTINUE
          END IF
        END IF
        jpt = 1
        IF (last<=kwi) THEN
          DO 80 j = last, lfirst, -1
            jpt = kwi - last + j
            cmb(jpt) = cmb(j)
80        CONTINUE
          DO 90 j = 1, jpt - 1
            cmb(j) = ' '
90        CONTINUE
        ELSE
          DO 100 j = lfirst, last
            jpt = kwi - last + j
            cmb(jpt) = cmb(j)
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
        CALL fmout(ma,cmb,kwi)
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
        CALL fmout(ma,cmb,kwi)
      ELSE
        GO TO 110
      END IF

      jform1 = jf1sav
      jform2 = jf2sav
      kflag = ksave
      RETURN

!             Error condition.

110   kflag = -8
      DO 120 j = 1, kwi
        cmb(j) = '*'
120   CONTINUE
      jform1 = jf1sav
      jform2 = jf2sav
      kflag = ksave
      RETURN
90000 FORMAT ('(I',I5,')')
    END SUBROUTINE zmfpcm
    SUBROUTINE zmfprt(form1,form2,ma)

!  Print MA in base 10 using FORM1 format for the real part and
!  FORM2 format for the imaginary part.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (*) :: form1, form2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz)
! ..
! .. Local Scalars ..
      INTEGER :: j, k, kwidim, kwidre, last, lsign
      CHARACTER (20) :: form
! ..
! .. External Subroutines ..
      EXTERNAL fmeq, zmfpcm
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, jformz, jprntz, kaccsw, &
        kdebug, keswch, kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, &
        ncall, ndg2mx, ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: m01(0:lunpck), m02(0:lunpck), m03(0:lunpck), m04(0:lunpck), &
        m05(0:lunpck), m06(0:lunpck), mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff), cmbufz(lmbufz)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fm1/m01, m02, m03, m04, m05, m06
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
      COMMON /zmbuff/cmbufz
      COMMON /zmuser/jformz, jprntz
! ..
      ncall = ncall + 1
      namest(ncall) = 'ZMFPRT'

      CALL zmfpcm(form1,ma,kwidre,cmbufz)
      CALL fmeq(ma(kptimu),m02)
      IF (m02(2)>=0) THEN
        lsign = 1
      ELSE
        lsign = -1
        m02(2) = -m02(2)
      END IF
      CALL zmfpcm(form2,m02,kwidim,cmbuff)

      cmbufz(kwidre+1) = ' '
      IF (lsign==1) THEN
        cmbufz(kwidre+2) = '+'
      ELSE
        cmbufz(kwidre+2) = '-'
      END IF
      cmbufz(kwidre+3) = ' '
      DO 10 j = 1, kwidim
        cmbufz(kwidre+3+j) = cmbuff(j)
10    CONTINUE
      cmbufz(kwidre+4+kwidim) = ' '
      cmbufz(kwidre+5+kwidim) = 'i'
      IF (jformz==2) cmbufz(kwidre+5+kwidim) = 'I'
      last = kwidre + kwidim + 5

      IF (m02(1)==mexpov .OR. m02(1)==mexpun) THEN
        DO 20 j = kwidre + 3, last
          IF (cmbufz(j)=='O' .OR. cmbufz(j)=='U') THEN
            cmbufz(j-2) = ' '
            GO TO 30
          END IF
20      CONTINUE
      END IF

30    WRITE (form,90000) kswide - 7
      WRITE (kw,form) (cmbufz(k),k=1,last)
      ncall = ncall - 1
      RETURN
90000 FORMAT (' (6X,',I3,'A1) ')
    END SUBROUTINE zmfprt
    SUBROUTINE zmi2m(integ,ma)

!  MA = INTEG

!  The real part of MA is set to the one word integer value INTEG.
!  The imaginary part is set to zero.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: integ
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz)
! ..
! .. External Subroutines ..
      EXTERNAL fmi2m, zmntr, zmntri
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
      namest(ncall) = 'ZMI2M '
      IF (ntrace/=0) CALL zmntri(2,integ,1)

      CALL fmi2m(integ,ma)
      CALL fmi2m(0,ma(kptimu))

      IF (ntrace/=0) CALL zmntr(1,ma,ma,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmi2m
    SUBROUTINE zm2i2m(integ1,integ2,ma)

!  MA = INTEG1 + INTEG2 i

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: integ1, integ2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz)
! ..
! .. External Subroutines ..
      EXTERNAL fmi2m, zmntr, zmntri
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
      namest(ncall) = 'ZM2I2M'
      IF (ntrace/=0) THEN
        CALL zmntri(2,integ1,1)
        CALL zmntri(2,integ2,0)
      END IF

      CALL fmi2m(integ1,ma)
      CALL fmi2m(integ2,ma(kptimu))

      IF (ntrace/=0) CALL zmntr(1,ma,ma,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zm2i2m
    SUBROUTINE zmimag(ma,mbfm)

!  MBFM = IMAG(MA)        imaginary part of MA

!  MA is a complex ZM number, MBFM is a real FM number.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mbfm(0:lunpck)
! ..
! .. External Subroutines ..
      EXTERNAL fmeq, fmntr, zmntr
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
      kflag = 0
      ncall = ncall + 1
      namest(ncall) = 'ZMIMAG'
      IF (ntrace/=0) CALL zmntr(2,ma,ma,1)

      CALL fmeq(ma(kptimu),mbfm)

      IF (ntrace/=0) CALL fmntr(1,mbfm,mbfm,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmimag
    SUBROUTINE zminp(line,ma,la,lb)

!  Convert an A1 character string to floating point multiple precision
!  complex format.

!  LINE is an A1 character array of length LB to be converted
!       to ZM format and returned in MA.
!  LA is a pointer telling the routine where in the array to begin
!     the conversion.  This allows more than one number to be stored
!     in an array and converted in place.
!  LB is a pointer to the last character of the field for that number.

!  The input numbers may be in integer or any real format.
!  In exponential format the 'E' may also be 'D', 'Q', or 'M'.

!  The following are all valid input strings:

!  1.23 + 4.56 I
!  1.23 + 4.56*I
!  2 + i
!  -i
!  1.23
!  4.56i
!  ( 1.23 , 4.56 )

!  So that ZMINP will convert any output from ZMOUT, LINE is tested
!  to see if the input contains any of the special symbols +OVERFLOW,
!  -OVERFLOW, +UNDERFLOW, -UNDERFLOW, or UNKNOWN.
!  For user input the abbreviations OVFL, UNFL, UNKN may be used.

      IMPLICIT NONE

!             Scratch array usage during ZMINP:   M01 - M05

!  Simulate a finite-state automaton to scan the input line
!  and build the number.  States 2-8 refer to the real part,
!  states 10-16 refer to the imaginary part.
!  States of the machine:

!   1.  Initial entry to the subroutine
!   2.  Sign of the number
!   3.  Scanning digits before a decimal point
!   4.  Decimal point
!   5.  Scanning digits after a decimal point
!   6.  E, D, Q, or M - precision indicator before the exponent
!   7.  Sign of the exponent
!   8.  Scanning exponent
!   9.  Comma between the real and imaginary part
!  10.  Sign of the number
!  11.  Scanning digits before a decimal point
!  12.  Decimal point
!  13.  Scanning digits after a decimal point
!  14.  E, D, Q, or M - precision indicator before the exponent
!  15.  Sign of the exponent
!  16.  Scanning exponent
!  17.  Syntax error

!  Character types recognized by the machine:

!  1.  Sign (+,-)
!  2.  Numeral (0,1,...,9)
!  3.  Decimal point (.)
!  4.  Precision indicator (E,D,Q,M)
!  5.  Illegal character for number
!  6.  Comma (,)
!  7.  Character to be ignored   ' '    '('    ')'    '*'

!  All blanks are ignored.  The analysis of the number proceeds as
!  follows:  If the simulated machine is in state JSTATE and a character
!  of type JTYPE is encountered the new state of the machine is given by
!  JTRANS(JSTATE,JTYPE).

!  State  1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16

! .. Intrinsic Functions ..
      INTRINSIC ichar, max, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: la, lb
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz)
      CHARACTER (1) :: line(lb)
! ..
! .. Local Scalars ..
      INTEGER :: j, jstate, k, kasave, kdigfl, kflag1, kiflag, kpt, krsave, &
        ksign, kstart, kstop, kstopi, kstopr, kstrti, kstrtr, ktype, kval, &
        ndsave, ntrsav
! ..
! .. Local Arrays ..
      INTEGER :: jtrans(16,4)
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmim, fminp, zmwarn
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
! .. Data Statements ..
      DATA jtrans/2, 17, 10, 10, 10, 7, 17, 10, 10, 17, 17, 17, 17, 15, 17, &
        17, 3, 3, 3, 5, 5, 8, 8, 8, 11, 11, 11, 13, 13, 16, 16, 16, 4, 4, 4, &
        17, 17, 17, 17, 17, 12, 12, 12, 17, 17, 17, 17, 17, 6, 6, 6, 6, 6, 8, &
        17, 17, 14, 14, 14, 14, 14, 16, 17, 17/
! ..
      IF (mblogs/=mbase) CALL fmcons
      ncall = ncall + 1
      namest(ncall) = 'ZMINP '
      ndsave = ndig
      kasave = kaccsw
      krsave = kround
      kround = 1
      kflag = 0

!             Since arithmetic tracing is not usually desired during
!             I/O conversion, disable tracing during this routine.

      ntrsav = ntrace
      ntrace = 0

!             Increase the working precision.

      IF (ncall<=2) THEN
        k = ngrd52
        ndig = max(ndig+k,2)
        IF (ndig>ndg2mx) THEN
          kflag = -9
          CALL zmwarn
          ma(0) = nint(ndig*alogm2)
          ma(1) = munkno
          ma(2) = 1
          ma(kptimu) = nint(ndig*alogm2)
          ma(kptimu+1) = munkno
          ma(kptimu+2) = 1
          DO 10 j = 2, ndsave
            ma(j+1) = 0
            ma(kptimu+j+1) = 0
10        CONTINUE
          GO TO 30
        END IF
      END IF
      kstart = la
      kstop = lb
      jstate = 1
      kstrtr = 0
      kstopr = 0
      kstrti = 0
      kstopi = 0
      kdigfl = 0
      kiflag = 0
      ksign = 1

!             Scan the number.

      DO 20 j = kstart, kstop
        IF (line(j)==' ' .OR. line(j)=='(' .OR. line(j)==')' .OR. &
          line(j)=='*') GO TO 20
        IF (line(j)=='I' .OR. line(j)=='i') THEN
          kiflag = 1
          IF (kstrti==0) THEN
            kstrti = kstrtr
            kstopi = kstopr
            kstrtr = 0
            kstopr = 0
          END IF
          GO TO 20
        END IF

        kpt = ichar(line(j))
        IF (kpt<lhash1 .OR. kpt>lhash2) THEN
          WRITE (kw,90000) line(j), kpt, lhash1, lhash2
          ktype = 5
          kval = 0
        ELSE
          ktype = khasht(kpt)
          kval = khashv(kpt)
        END IF
        IF (ktype==2 .OR. ktype==5) kdigfl = 1
        IF (line(j)==',') THEN
          IF (jstate<9) THEN
            jstate = 9
          ELSE
            GO TO 40
          END IF
        ELSE
          IF (ktype>=5) ktype = 2
          IF (jstate<17) jstate = jtrans(jstate,ktype)
        END IF
        IF (jstate==9 .OR. jstate==10) kdigfl = 0
        IF (jstate==2 .OR. jstate==10) ksign = kval

        IF (jstate>=2 .AND. jstate<=8) THEN
          IF (kstrtr==0) kstrtr = j
          kstopr = j
        END IF
        IF (jstate>=10 .AND. jstate<=16) THEN
          IF (kstrti==0) kstrti = j
          kstopi = j
        END IF

20    CONTINUE

!             Form the number and return.

      IF (kstrtr>0) THEN
        CALL fminp(line,ma,kstrtr,kstopr)
      ELSE
        CALL fmim(0,ma)
      END IF
      kflag1 = kflag

      IF (kstrti>0) THEN
        IF (kiflag==1 .AND. kdigfl==0) THEN
          CALL fmim(ksign,ma(kptimu))
        ELSE
          CALL fminp(line,ma(kptimu),kstrti,kstopi)
        END IF
      ELSE IF (kiflag==1) THEN
        CALL fmim(1,ma(kptimu))
      ELSE
        CALL fmim(0,ma(kptimu))
      END IF

      IF (kflag1/=0 .OR. kflag/=0 .OR. jstate==17) GO TO 40

30    ndig = ndsave
      kaccsw = kasave
      ntrace = ntrsav
      kround = krsave
      IF (kflag==1) kflag = 0
      ma(0) = nint(ndig*alogm2)
      ma(kptimu) = ma(0)
      ncall = ncall - 1
      RETURN

!             Error in converting the number.

40    kflag = -7
      CALL zmwarn
      ma(0) = nint(ndig*alogm2)
      ma(1) = munkno
      ma(2) = 1
      ma(kptimu) = nint(ndig*alogm2)
      ma(kptimu+1) = munkno
      ma(kptimu+2) = 1
      DO 50 j = 2, ndsave
        ma(j+1) = 0
        ma(kptimu+j+1) = 0
50    CONTINUE
      GO TO 30
90000 FORMAT (/' Error in input conversion.'/ &
        ' ICHAR function was out of range for the current', &
        ' dimensions.'/' ICHAR(''',A,''') gave the value ',I12, &
        ', which is outside the currently'/' dimensioned',' bounds of (',I5, &
        ':',I5,') for variables KHASHT ','and KHASHV.'/ &
        ' Re-define the two parameters ', &
        'LHASH1 and LHASH2 so the dimensions will'/' contain', &
        ' all possible output values from ICHAR.'//)
    END SUBROUTINE zminp
    SUBROUTINE zmint(ma,mb)

!  MB = INT(MA)

!  The integer parts of both real and imaginary values are returned.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. External Subroutines ..
      EXTERNAL fmint, zmntr
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
      namest(ncall) = 'ZMINT '
      IF (ntrace/=0) CALL zmntr(2,ma,ma,1)

      CALL fmint(ma,mb)
      CALL fmint(ma(kptimu),mb(kptimu))

      IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmint
    SUBROUTINE zmipwr(ma,ival,mb)

!  MB = MA ** IVAL

!  Raise a ZM number to an integer power.
!  The binary multiplication method used requires an average of
!  1.5 * LOG2(IVAL) multiplications.

      IMPLICIT NONE

!             Scratch array usage during ZMIPWR:  M01 - M03, MZ01 - MZ02

! .. Intrinsic Functions ..
      INTRINSIC abs, int, log, max, min, mod, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ma2, maccmb, maiz, marz, mxsave
      REAL :: xval
      INTEGER :: i2n, j, k, kasave, kovun, kwrnsv, lvlsav, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmcons, fmim, fmipwr, fmntri, zmdiv, zmeq, zmeq2, zmexit, &
        zmi2m, zmmpy, zmntr, zmsqr, zmwarn
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
        mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), mz03(0:lunpkz), &
        mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      ncall = ncall + 1
      namest(ncall) = 'ZMIPWR'
      ndsave = ndig
      IF (ntrace/=0) THEN
        CALL zmntr(2,ma,ma,1)
        CALL fmntri(2,ival,0)
      END IF
      kovun = 0
      marz = ma(0)
      maiz = ma(kptimu)
      IF (ma(1)==mexpov .OR. ma(1)==mexpun .OR. ma(kptimu+1)==mexpov .OR. &
        ma(kptimu+1)==mexpun) kovun = 1

      IF (mblogs/=mbase) CALL fmcons
      kflag = 0
      kasave = kaccsw
      mxsave = mxexp
      mxexp = mxexp2

!             Check for special cases.

      IF (ma(1)==munkno .OR. ma(kptimu+1)==munkno .OR. (ival<=0 .AND. ma( &
          2)==0 .AND. ma(kptimu+2)==0)) THEN
        ma2 = ma(2)
        mb(0) = nint(ndig*alogm2)
        mb(1) = munkno
        mb(2) = 1
        mb(kptimu) = nint(ndig*alogm2)
        mb(kptimu+1) = munkno
        mb(kptimu+2) = 1
        DO 10 j = 2, ndsave
          mb(j+1) = 0
          mb(kptimu+j+1) = 0
10      CONTINUE
        kflag = -4
        IF (ival<=0 .AND. ma2==0) CALL zmwarn
        IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
        ncall = ncall - 1
        mxexp = mxsave
        RETURN
      END IF

      IF (ival==0) THEN
        CALL zmi2m(1,mb)
        IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
        ncall = ncall - 1
        mxexp = mxsave
        RETURN
      END IF

      IF (abs(ival)==1) THEN
        kwrnsv = kwarn
        kwarn = 0
        IF (ival==1) THEN
          CALL zmeq(ma,mb)
        ELSE
          k = int((5.0D0*dlogtn)/dlogmb+2.0D0)
          ndig = min(max(ndig+k,2),ndg2mx)
          CALL zmi2m(1,mz02)
          CALL zmeq2(ma,ma,ndsave,ndig,1)
          CALL zmdiv(mz02,ma,mb)
          CALL zmeq2(mb,mb,ndig,ndsave,1)
          ndig = ndsave
        END IF
        IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
        ncall = ncall - 1
        kwarn = kwrnsv
        mxexp = mxsave
        RETURN
      END IF

      IF (ma(2)==0 .AND. ma(kptimu+2)==0) THEN
        CALL zmi2m(0,mb)
        IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
        ncall = ncall - 1
        mxexp = mxsave
        RETURN
      END IF

      IF (ma(kptimu+2)==0) THEN
        ncall = ncall - 1
        lvlsav = lvltrc
        lvltrc = lvltrc - 1
        CALL fmipwr(ma,ival,mb)
        CALL fmim(0,mb(kptimu))
        ncall = ncall + 1
        lvltrc = lvlsav
        IF (ntrace/=0) THEN
          namest(ncall) = 'ZMIPWR'
          CALL zmntr(1,mb,mb,1)
        END IF
        ncall = ncall - 1
        mxexp = mxsave
        RETURN
      END IF

      IF (ma(2)==0) THEN
        ncall = ncall - 1
        lvlsav = lvltrc
        lvltrc = lvltrc - 1
        IF (ival>=0) THEN
          i2n = mod(ival,4)
        ELSE
          i2n = mod(4-mod(abs(ival),4),4)
        END IF
        IF (i2n==0) THEN
          CALL fmipwr(ma(kptimu),ival,mb)
          CALL fmim(0,mb(kptimu))
        ELSE IF (i2n==1) THEN
          CALL fmipwr(ma(kptimu),ival,mb(kptimu))
          CALL fmim(0,mb)
        ELSE IF (i2n==2) THEN
          CALL fmipwr(ma(kptimu),ival,mb)
          CALL fmim(0,mb(kptimu))
          IF (mb(1)/=munkno) mb(2) = -mb(2)
        ELSE IF (i2n==3) THEN
          CALL fmipwr(ma(kptimu),ival,mb(kptimu))
          CALL fmim(0,mb)
          IF (mb(kptimu+1)/=munkno) mb(kptimu+2) = -mb(kptimu+2)
        END IF
        ncall = ncall + 1
        lvltrc = lvlsav
        IF (ntrace/=0) THEN
          namest(ncall) = 'ZMIPWR'
          CALL zmntr(1,mb,mb,1)
        END IF
        ncall = ncall - 1
        mxexp = mxsave
        RETURN
      END IF

!             Increase the working precision.

      IF (ncall==1) THEN
        xval = abs(ival) + 1
        k = int((5.0*real(dlogtn)+1.5*log(xval))/alogmb+2.0)
        ndig = max(ndig+k,2)
      ELSE
        xval = abs(ival) + 1
        k = int(log(xval)/alogmb+1.0)
        ndig = ndig + k
      END IF
      IF (ndig>ndg2mx) THEN
        kflag = -9
        CALL zmwarn
        mb(1) = munkno
        mb(2) = 1
        mb(kptimu+1) = munkno
        mb(kptimu+2) = 1
        DO 20 j = 2, ndsave
          mb(j+1) = 0
          mb(kptimu+j+1) = 0
20      CONTINUE
        ndig = ndsave
        mb(0) = nint(ndig*alogm2)
        mb(kptimu) = nint(ndig*alogm2)
        ndig = ndsave
        IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
        mxexp = mxsave
        kaccsw = kasave
        ncall = ncall - 1
        RETURN
      END IF

!             Initialize.

      kwrnsv = kwarn
      kwarn = 0
      k = abs(ival)

      CALL zmeq2(ma,mz02,ndsave,ndig,0)

      IF (mod(k,2)==0) THEN
        CALL zmi2m(1,mb)
      ELSE
        CALL zmeq(mz02,mb)
      END IF

!             This is the multiplication loop.

30    k = k/2
      CALL zmsqr(mz02,mz02)
      IF (mod(k,2)==1) CALL zmmpy(mz02,mb,mb)
      IF (k>1) GO TO 30

!             Invert if the exponent is negative.

      IF (ival<0) THEN
        CALL zmi2m(1,mz02)
        CALL zmdiv(mz02,mb,mb)
      END IF
      kwarn = kwrnsv

!             Round the result and return.

      maccmb = mb(0)
      ma(0) = marz
      mb(0) = min(maccmb,marz,maiz)
      maccmb = mb(kptimu)
      ma(kptimu) = maiz
      mb(kptimu) = min(maccmb,marz,maiz)
      CALL zmexit(mb,mb,ndsave,mxsave,kasave,kovun,1)
      RETURN
    END SUBROUTINE zmipwr
    SUBROUTINE zmlg10(ma,mb)

!  MB = LOG10(MA).

      IMPLICIT NONE

!             Scratch array usage during ZMLG10:  M01 - M05, MZ01 - MZ02

! .. Intrinsic Functions ..
      INTRINSIC min
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccmb, maiz, marz, mxsave
      INTEGER :: kasave, kovun, kreslt, krsave, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmdivd, fmlni, zmentr, zmeq2, zmexit, zmln
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
        mpisav(0:lunpck), mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      CALL zmentr('ZMLG10',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
      IF (kreslt/=0) RETURN
      marz = ma(0)
      maiz = ma(kptimu)
      kaccsw = 0
      krsave = krad
      krad = 1

      CALL zmeq2(ma,ma,ndsave,ndig,1)
      CALL zmln(ma,mz02)
      CALL fmlni(10,m03)
      CALL fmdivd(mz02,mz02(kptimu),m03,mz01,mz01(kptimu))

      maccmb = mz01(0)
      ma(0) = marz
      mz01(0) = min(maccmb,marz,maiz)
      maccmb = mz01(kptimu)
      ma(kptimu) = maiz
      mz01(kptimu) = min(maccmb,marz,maiz)
      CALL zmexit(mz01,mb,ndsave,mxsave,kasave,kovun,0)
      krad = krsave
      RETURN
    END SUBROUTINE zmlg10
    SUBROUTINE zmln(ma,mb)

!  MB = LN(MA).

      IMPLICIT NONE

!             Scratch array usage during ZMLN:   M01 - M05, MZ01

! .. Intrinsic Functions ..
      INTRINSIC int, min, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccmb, maiz, marz, mxsave
      INTEGER :: j, kasave, kf1, kovun, kreslt, krsave, ndsave
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL fmabs, fmadd, fmatn2, fmdiv, fmdivi, fmeq, fmi2m, fmln, fmmpy, &
        fmpi, fmsqr, fmsub, zmabs, zmentr, zmeq2, zmexit, zmntr, zmrslt, &
        zmwarn
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
        mpisav(0:lunpck), mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      CALL zmentr('ZMLN  ',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
      IF (kreslt/=0) RETURN
      marz = ma(0)
      maiz = ma(kptimu)
      kaccsw = 0
      krsave = krad
      krad = 1

      CALL zmeq2(ma,ma,ndsave,ndig,1)

!             Check for special cases.

      IF (ma(2)==0 .AND. ma(kptimu+2)==0) THEN
        kflag = -4
        mz01(1) = munkno
        mz01(2) = 1
        mz01(kptimu+1) = munkno
        mz01(kptimu+2) = 1
        DO 10 j = 2, ndsave
          mz01(j+1) = 0
          mz01(kptimu+j+1) = 0
10      CONTINUE
        mz01(0) = nint(ndig*alogm2)
        mz01(kptimu) = nint(ndig*alogm2)
        GO TO 20
      ELSE IF (ma(kptimu+2)==0) THEN
        IF (ma(2)<0) THEN
          CALL fmeq(ma,mz01)
          IF (mz01(1)/=munkno) mz01(2) = -mz01(2)
          CALL fmln(mz01,mz01)
          CALL fmpi(mz01(kptimu))
        ELSE
          CALL fmln(ma,mz01)
          CALL fmi2m(0,mz01(kptimu))
        END IF
        GO TO 20
      ELSE IF (ma(2)==0) THEN
        IF (ma(kptimu+2)<0) THEN
          CALL fmeq(ma(kptimu),mz01)
          IF (mz01(1)/=munkno) mz01(2) = -mz01(2)
          CALL fmln(mz01,mz01)
          CALL fmpi(mz01(kptimu))
          CALL fmdivi(mz01(kptimu),-2,mz01(kptimu))
        ELSE
          CALL fmln(ma(kptimu),mz01)
          CALL fmpi(mz01(kptimu))
          CALL fmdivi(mz01(kptimu),2,mz01(kptimu))
        END IF
        GO TO 20
      END IF

!             Ln(a + b i) = Ln(Abs(a + b i)) + Arg(a + b i) i.

      CALL fmabs(ma,m03)
      CALL fmabs(ma(kptimu),m04)

!             Check for cancellation in Ln(x).

      CALL fmi2m(1,m05)
      kf1 = 0
      IF (fmcomp(m03,'EQ',m05) .AND. m04(1)<=(-ndig)) kf1 = 1
      IF (fmcomp(m04,'EQ',m05) .AND. m03(1)<=(-ndig)) kf1 = 1

      IF (fmcomp(m03,'GE',m04)) THEN
        CALL fmsub(ma,m05,m03)
        CALL fmadd(ma,m05,m04)
        CALL fmmpy(m03,m04,m03)
        CALL fmsqr(ma(kptimu),m04)
        CALL fmadd(m03,m04,m04)
      ELSE
        CALL fmsub(ma(kptimu),m05,m03)
        CALL fmadd(ma(kptimu),m05,m04)
        CALL fmmpy(m03,m04,m03)
        CALL fmsqr(ma,m04)
        CALL fmadd(m03,m04,m04)
      END IF
      CALL zmabs(ma,mz01)
      CALL fmadd(mz01,m05,m03)
      CALL fmdiv(m04,m03,m03)
      IF (kf1==1) THEN
        CALL fmeq(m03,mz01)
        CALL fmatn2(ma(kptimu),ma,mz01(kptimu))
        GO TO 20
      ELSE IF (m03(1)<0) THEN
        ndig = ndig - int(m03(1))
        IF (ndig>ndg2mx) THEN
          namest(ncall) = 'ZMLN  '
          kflag = -9
          CALL zmwarn
          kreslt = 12
          ndig = ndsave
          CALL zmrslt(mb,kreslt)
          IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
          ncall = ncall - 1
          mxexp = mxsave
          kaccsw = kasave
          RETURN
        END IF
        CALL zmeq2(ma,ma,ndsave,ndig,1)
        CALL zmabs(ma,mz01)
      END IF

      CALL fmln(mz01,mz01)
      CALL fmatn2(ma(kptimu),ma,mz01(kptimu))

20    maccmb = mz01(0)
      ma(0) = marz
      mz01(0) = min(maccmb,marz,maiz)
      maccmb = mz01(kptimu)
      ma(kptimu) = maiz
      mz01(kptimu) = min(maccmb,marz,maiz)
      CALL zmexit(mz01,mb,ndsave,mxsave,kasave,kovun,0)
      krad = krsave
      RETURN
    END SUBROUTINE zmln
    SUBROUTINE zmm2i(ma,integ)

!  INTEG = MA

!  INTEG is set to the integer value of the real part of MA

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: integ
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz)
! ..
! .. External Subroutines ..
      EXTERNAL fmm2i, zmntr, zmntri
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
      namest(ncall) = 'ZMM2I '
      IF (ntrace/=0) CALL zmntr(2,ma,ma,1)

      CALL fmm2i(ma,integ)

      IF (ntrace/=0) CALL zmntri(1,integ,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmm2i
    SUBROUTINE zmm2z(ma,zval)

!  ZVAL = MA

!  Complex variable ZVAL is set to MA.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC cmplx
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      COMPLEX :: zval
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL :: di, dr
! ..
! .. External Subroutines ..
      EXTERNAL fmm2sp, zmntr, zmntrz
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
      namest(ncall) = 'ZMM2Z '
      IF (ntrace/=0) CALL zmntr(2,ma,ma,1)

      CALL fmm2sp(ma,dr)
      CALL fmm2sp(ma(kptimu),di)
      zval = cmplx(dr,di)

      IF (ntrace/=0) CALL zmntrz(1,zval,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmm2z
    SUBROUTINE zmmpy(ma,mb,mc)

!  MC = MA * MB

      IMPLICIT NONE

!             Scratch array usage during ZMMPY:   M01 - M03, MZ01

! .. Intrinsic Functions ..
      INTRINSIC abs, int, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz), mc(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccmb, maiz, marz, mbiz, mbrz, mxsave, mz11sv, mzero
      INTEGER :: iextra, j, kasave, kmethd, kovun, kreslt, kwrnsv, ndgsv2, &
        ndsave, ngoal, ntrsav
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmmpy, fmmpyd, fmsub, zmentr, zmeq2, zmntr, zmrslt, &
        zmwarn
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
        mpisav(0:lunpck), mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      IF (abs(ma(1))>mexpab .OR. abs(ma(kptimu+1))>mexpab .OR. abs(mb(1))> &
          mexpab .OR. abs(mb(kptimu+1))>mexpab .OR. kdebug>=1) THEN
        CALL zmentr('ZMMPY ',ma,mb,2,mc,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        IF (ntrace/=0) THEN
          namest(ncall) = 'ZMMPY '
          CALL zmntr(2,ma,mb,2)
        END IF
        ndsave = ndig
        IF (ncall==1) THEN
          ndig = max(ndig+ngrd52,2)
          IF (ndig>ndg2mx) THEN
            namest(ncall) = 'ZMMPY '
            kflag = -9
            CALL zmwarn
            kreslt = 12
            ndig = ndsave
            CALL zmrslt(mc,kreslt)
            IF (ntrace/=0) CALL zmntr(1,mc,mc,1)
            ncall = ncall - 1
            RETURN
          END IF
          IF (mbase>=100*abs(ma(2)) .OR. mbase>=100*abs(ma(kptimu+2))) THEN
            ndig = min(ndig+1,ndg2mx)
          ELSE IF (mbase>=100*abs(mb(2)) .OR. mbase>=100*abs(mb(kptimu+ &
              2))) THEN
            ndig = min(ndig+1,ndg2mx)
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 1
        mxsave = mxexp
        mxexp = mxexp2
        kovun = 0
      END IF

      marz = ma(0)
      mbrz = mb(0)
      maiz = ma(kptimu)
      mbiz = mb(kptimu)
      mz11sv = -munkno
      mzero = 0
      ntrsav = ntrace
      ntrace = 0
      kwrnsv = kwarn
      kwarn = 0

10    DO 20 j = ndsave + 2, ndig + 1
        ma(j) = mzero
        mb(j) = mzero
        ma(kptimu+j) = mzero
        mb(kptimu+j) = mzero
20    CONTINUE
      IF (ncall==1) THEN
        ma(0) = nint(ndig*alogm2)
        mb(0) = ma(0)
        ma(kptimu) = ma(0)
        mb(kptimu) = ma(0)
      END IF

!             Check for special cases.

      kmethd = 1
      IF (ndig>=35) kmethd = 2

      IF (mb(kptimu+2)==0) THEN
        CALL fmmpyd(mb,ma,ma(kptimu),mz01,mz01(kptimu))
      ELSE IF (mb(2)==0) THEN
        CALL fmmpyd(mb(kptimu),ma(kptimu),ma,mz01,mz01(kptimu))
        IF (mz01(1)/=munkno) mz01(2) = -mz01(2)
      ELSE IF (ma(kptimu+2)==0) THEN
        CALL fmmpyd(ma,mb,mb(kptimu),mz01,mz01(kptimu))
      ELSE IF (ma(2)==0) THEN
        CALL fmmpyd(ma(kptimu),mb(kptimu),mb,mz01,mz01(kptimu))
        IF (mz01(1)/=munkno) mz01(2) = -mz01(2)
      ELSE IF (kmethd==1) THEN

!             Method 1 for  ( a + b i ) * ( c + d i )

!             result = a*c - b*d + ( a*d + b*c ) i

        kaccsw = 0
        CALL fmmpyd(ma,mb,mb(kptimu),mz01,mz01(kptimu))
        CALL fmmpyd(ma(kptimu),mb(kptimu),mb,m01,m02)
        IF (mz01(2)*m01(2)<0) THEN
          kaccsw = 0
        ELSE
          kaccsw = 1
        END IF
        CALL fmsub(mz01,m01,mz01)
        IF (mz01(kptimu+2)*m02(2)<0) THEN
          kaccsw = 1
        ELSE
          kaccsw = 0
        END IF
        CALL fmadd(mz01(kptimu),m02,mz01(kptimu))
        kaccsw = 1
      ELSE

!             Method 2 for  ( a + b i ) * ( c + d i )

!             P = ( a + b )*( c + d )
!             result = a*c - b*d + ( P - a*c - b*d ) i

        CALL fmadd(ma,ma(kptimu),m01)
        CALL fmadd(mb,mb(kptimu),m02)
        CALL fmmpy(m01,m02,m01)

        CALL fmmpy(ma,mb,m02)
        CALL fmmpy(ma(kptimu),mb(kptimu),m03)

        CALL fmsub(m02,m03,mz01)
        CALL fmsub(m01,m02,mz01(kptimu))
        CALL fmsub(mz01(kptimu),m03,mz01(kptimu))
      END IF

!             Check for too much cancellation.

      IF (ncall<=1) THEN
        ngoal = int(real(ndsave)*alogm2) + 7
      ELSE
        ngoal = int(-mxexp2)
      END IF
      IF (mz01(0)<=ngoal .OR. mz01(kptimu)<=ngoal) THEN
        IF (mz11sv>-munkno .AND. mz01(0)>ngoal .AND. mz01(kptimu+2)==0) &
          GO TO 40
        IF (mz11sv>-munkno .AND. mz01(kptimu)>ngoal .AND. mz01(2)==0) GO TO 40
        iextra = int(real(max(ngoal-mz01(0),ngoal-mz01(kptimu)))/alogm2+23.03/ &
          alogmb) + 1
        ndig = ndig + iextra
        IF (ndig>ndg2mx) THEN
          namest(ncall) = 'ZMMPY '
          kflag = -9
          CALL zmwarn
          mz01(1) = munkno
          mz01(2) = 1
          mz01(kptimu+1) = munkno
          mz01(kptimu+2) = 1
          DO 30 j = 2, ndsave
            mz01(j+1) = 0
            mz01(kptimu+j+1) = 0
30        CONTINUE
          ndig = ndig - iextra
          mz01(0) = nint(ndig*alogm2)
          mz01(kptimu) = nint(ndig*alogm2)
          GO TO 40
        END IF
        mz11sv = mz01(1)
        GO TO 10
      END IF

40    mxexp = mxsave
      ntrace = ntrsav
      ndgsv2 = ndig
      ndig = ndsave
      kwarn = kwrnsv
      maccmb = mz01(0)
      ma(0) = marz
      mb(0) = mbrz
      mz01(0) = min(maccmb,marz,maiz,mbrz,mbiz)
      maccmb = mz01(kptimu)
      ma(kptimu) = maiz
      mb(kptimu) = mbiz
      mz01(kptimu) = min(maccmb,marz,maiz,mbrz,mbiz)
      CALL zmeq2(mz01,mc,ndgsv2,ndsave,0)
      IF (mc(1)>=mexpov .OR. mc(1)<=-mexpov .OR. mc(kptimu+1)>=mexpov .OR. &
          mc(kptimu+1)<=-mexpov) THEN
        IF (mc(1)==munkno .OR. mc(kptimu+1)==munkno) THEN
          kflag = -4
        ELSE IF (mc(1)==mexpov .OR. mc(kptimu+1)==mexpov) THEN
          kflag = -5
        ELSE IF (mc(1)==mexpun .OR. mc(kptimu+1)==mexpun) THEN
          kflag = -6
        END IF
        IF ((mc(1)==munkno) .OR. (mc(kptimu+1)==munkno) .OR. (mc(1)==mexpun &
            .AND. kovun==0) .OR. (mc(kptimu+1)==mexpun .AND. kovun==0) .OR. ( &
            mc(1)==mexpov .AND. kovun==0) .OR. (mc(kptimu+ &
            1)==mexpov .AND. kovun==0)) THEN
          namest(ncall) = 'ZMMPY '
          CALL zmwarn
        END IF
      END IF
      IF (ntrace/=0) CALL zmntr(1,mc,mc,1)
      kaccsw = kasave
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmmpy
    SUBROUTINE zmmpyi(ma,integ,mb)

!  MB = MA * INTEG        Multiply by one-word (real) integer.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: integ
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: mxsave
      INTEGER :: kasave, kovun, kreslt, kwrnsv, ndsave, ntrsav
! ..
! .. External Subroutines ..
      EXTERNAL fmmpyi, fmntri, zmentr, zmntr, zmwarn
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
      IF (abs(ma(1))>mexpab .OR. abs(ma(kptimu+1))>mexpab .OR. kdebug>=1) THEN
        CALL zmentr('ZMMPYI',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
        ndig = ndsave
        mxexp = mxsave
        kaccsw = kasave
        ntrsav = ntrace
      ELSE
        ncall = ncall + 1
        IF (ntrace/=0) THEN
          namest(ncall) = 'ZMMPYI'
          CALL zmntr(2,ma,ma,1)
          CALL fmntri(2,integ,0)
        END IF
        kovun = 0
      END IF

!             Force FMMPYI to use more guard digits for user calls.

      ncall = ncall - 1
      ntrsav = ntrace
      ntrace = 0
      kwrnsv = kwarn
      kwarn = 0

      CALL fmmpyi(ma,integ,mb)
      CALL fmmpyi(ma(kptimu),integ,mb(kptimu))

      ntrace = ntrsav
      kwarn = kwrnsv
      ncall = ncall + 1
      IF (ntrace/=0) namest(ncall) = 'ZMMPYI'
      IF (mb(1)==munkno .OR. mb(kptimu+1)==munkno) THEN
        kflag = -4
      ELSE IF (mb(1)==mexpov .OR. mb(kptimu+1)==mexpov) THEN
        kflag = -5
      ELSE IF (mb(1)==mexpun .OR. mb(kptimu+1)==mexpun) THEN
        kflag = -6
      END IF
      IF ((mb(1)==munkno) .OR. (mb(kptimu+1)==munkno) .OR. (mb(1)==mexpun &
          .AND. kovun==0) .OR. (mb(kptimu+1)==mexpun .AND. kovun==0) .OR. (mb( &
          1)==mexpov .AND. kovun==0) .OR. (mb(kptimu+ &
          1)==mexpov .AND. kovun==0)) THEN
        namest(ncall) = 'ZMMPYI'
        CALL zmwarn
      END IF
      IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmmpyi
    SUBROUTINE zmnint(ma,mb)

!  MB = NINT(MA)

!  The nearest integers to both real and imaginary parts are returned.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. External Subroutines ..
      EXTERNAL fmnint, zmntr
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
      namest(ncall) = 'ZMNINT'
      IF (ntrace/=0) CALL zmntr(2,ma,ma,1)

      CALL fmnint(ma,mb)
      CALL fmnint(ma(kptimu),mb(kptimu))

      IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmnint
    SUBROUTINE zmntr(ntr,ma,mb,narg)

!  Print ZM numbers in base 10 format using ZMOUT for conversion.
!  This is used for trace output from the ZM routines.

!  NTR =  1 if a result of an ZM call is to be printed.
!      =  2 to print input argument(s) to an ZM call.

!  MA  -  the ZM number to be printed.

!  MB  -  an optional second ZM number to be printed.

!  NARG - the number of arguments.  NARG = 1 if only MA is to be
!         printed, and NARG = 2 if both MA and MB are to be printed.

!  NTRACE and LVLTRC (in COMMON /FMUSER/) control trace printout.

!  NTRACE = 0        No printout except warnings and errors.

!  NTRACE = 1        The result of each call to one of the routines
!                    is printed in base 10, using ZMOUT.

!  NTRACE = -1       The result of each call to one of the routines
!                    is printed in internal base MBASE format.

!  NTRACE = 2        The input arguments and result of each call to one
!                    of the routines is printed in base 10, using ZMOUT.

!  NTRACE = -2       The input arguments and result of each call to one
!                    of the routines is printed in base MBASE format.

!  LVLTRC defines the call level to which the trace is done.  LVLTRC = 1
!         means only FM routines called directly by the user are traced,
!         LVLTRC = K prints traces for ZM or FM routines with call
!         levels up to and including level K.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, int
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: narg, ntr
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. Local Scalars ..
      CHARACTER (6) :: name
! ..
! .. External Subroutines ..
      EXTERNAL zmntrj, zmprnt
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
        CALL zmntrj(ma,ndig)
        IF (narg==2) CALL zmntrj(mb,ndig)
      END IF

!             Check for base 10 trace using ZMOUT.

      IF (ntrace>0) THEN
        CALL zmprnt(ma)

        IF (narg==2) THEN
          CALL zmprnt(mb)
        END IF
      END IF

      RETURN
90000 FORMAT (' Input to ',A6)
90010 FORMAT (' ',A6,15X,'Call level =',I2,5X,'MBASE =',I10,5X,'NDIG =',I6)
90020 FORMAT (' ',A6,6X,'Call level =',I2,4X,'MBASE =',I10,4X,'NDIG =',I6,4X, &
        'KFLAG =',I3)
    END SUBROUTINE zmntr
    SUBROUTINE zmntr2(ntr,mafm,mbfm,narg)

!  Print real FM numbers in base 10 format using FMOUT for conversion.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, int
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: narg, ntr
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: mafm(0:lunpck), mbfm(0:lunpck)
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
        CALL fmntrj(mafm,ndig)
        IF (narg==2) CALL fmntrj(mbfm,ndig)
      END IF

!             Check for base 10 trace using FMOUT.

      IF (ntrace>0) THEN
        CALL fmprnt(mafm)

        IF (narg==2) THEN
          CALL fmprnt(mbfm)
        END IF
      END IF

      RETURN
90000 FORMAT (' Input to ',A6)
90010 FORMAT (' ',A6,15X,'Call level =',I2,5X,'MBASE =',I10,5X,'NDIG =',I6)
90020 FORMAT (' ',A6,6X,'Call level =',I2,4X,'MBASE =',I10,4X,'NDIG =',I6,4X, &
        'KFLAG =',I3)
    END SUBROUTINE zmntr2
    SUBROUTINE zmntri(ntr,n,knam)

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
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
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
    END SUBROUTINE zmntri
    SUBROUTINE zmntrj(ma,nd)

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
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: nd
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz)
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
      WRITE (kw,form) (int(ma(j+kptimu)),j=1,n1)

      RETURN
90000 FORMAT (' (1X,I19,I',I2,',',I3,'I',I2,') ')
90010 FORMAT (' (1X,I19,I',I2,',',I3,'I',I2,'/(22X,',I3,'I',I2,')) ')
    END SUBROUTINE zmntrj
    SUBROUTINE zmntrz(ntr,x,knam)

!  Internal routine for trace output of complex variables.

!  NTR - 1 for output values
!        2 for input values

!  X   - Complex value to be printed if NX.EQ.1

!  KNAM - Positive if the routine name is to be printed.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs, aimag, dble, int
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      COMPLEX :: x
      INTEGER :: knam, ntr
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ximag, xreal
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

      xreal = dble(x)
      ximag = dble(aimag(x))
      IF (ximag>=0.0D0) THEN
        WRITE (kw,90030) xreal, ximag
      ELSE
        WRITE (kw,90040) xreal, abs(ximag)
      END IF

      RETURN
90000 FORMAT (' Input to ',A6)
90010 FORMAT (' ',A6,15X,'Call level =',I2,5X,'MBASE =',I10,5X,'NDIG =',I6)
90020 FORMAT (' ',A6,6X,'Call level =',I2,4X,'MBASE =',I10,4X,'NDIG =',I6,4X, &
        'KFLAG =',I3)
90030 FORMAT (1X,D30.20,' +',D30.20,' i')
90040 FORMAT (1X,D30.20,' -',D30.20,' i')
    END SUBROUTINE zmntrz
    SUBROUTINE zmout(ma,line,lb,last1,last2)

!  Convert a floating multiple precision number to a character array
!  for output.

!  MA    is an ZM number to be converted to an A1 character
!        array in base 10 format
!  LINE  is the CHARACTER*1 array in which the result is returned.
!  LB    is the length of LINE.
!  LAST1 is the position of the last nonblank character of the
!        real part of the number in LINE.
!  LAST2 is the position of the last nonblank character of the
!        imaginary part of the number in LINE.

!  JFORM1 and JFORM2 determine the format of the two FM numbers
!  making up the complex value MA.  See FMOUT for details.

!  JFORMZ determines the format of the real and imaginary parts.

!  JFORMZ = 1  normal setting :       1.23 - 4.56 i
!         = 2  use capital I  :       1.23 - 4.56 I
!         = 3  parenthesis format   ( 1.23 , -4.56 )

!  LINE should be dimensioned at least 4*(LOG10(MBASE)*NDIG + 15) on a
!  32-bit machine to allow for up to 10 digit exponents.  Replace
!  15 by 20 if 48-bit integers are used, 25 for 64-bit integers, etc.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC int, log10, max, min, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: last1, last2, lb
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz)
      CHARACTER (1) :: line(lb)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maim2
      INTEGER :: j, kpt, lb2, nd, nexp
! ..
! .. External Subroutines ..
      EXTERNAL fmout
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, jformz, jprntz, kaccsw, &
        kdebug, keswch, kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, &
        ncall, ndg2mx, ndig, ntrace
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
      COMMON /zmuser/jformz, jprntz
! ..
      ncall = ncall + 1
      namest(ncall) = 'ZMOUT '
      DO 10 j = 1, lb
        line(j) = ' '
10    CONTINUE
      nd = int(real(ndig)*log10(real(mbase))) + 1
      IF (nd<2) nd = 2
      nexp = int(2.0*log10(real(mxbase))) + 6
      kpt = 1
      IF (jformz==3) kpt = 3
      lb2 = max(jform2+nexp,nd+nexp)
      lb2 = min(lb+1-kpt,lb2)
      CALL fmout(ma,line(kpt),lb2)

      IF (jformz==3) line(1) = '('
      last1 = 1
      DO 20 j = lb2, 1, -1
        IF (line(j)/=' ') THEN
          last1 = j
          GO TO 30
        END IF
20    CONTINUE

30    maim2 = ma(kptimu+2)
      line(last1+1) = ' '
      IF (jformz==3) THEN
        line(last1+2) = ','
      ELSE
        IF (maim2<0) THEN
          ma(kptimu+2) = -ma(kptimu+2)
          line(last1+2) = '-'
        ELSE
          line(last1+2) = '+'
        END IF
      END IF

      kpt = last1 + 3
      lb2 = max(jform2+nexp,nd+nexp)
      lb2 = min(lb+1-kpt,lb2+2)
      CALL fmout(ma(kptimu),line(kpt),lb2)
      last1 = kpt
      DO 40 j = lb2 + kpt - 1, kpt, -1
        IF (line(j)/=' ') THEN
          last2 = j
          GO TO 50
        END IF
40    CONTINUE

50    last2 = last2 + 2
      line(last2) = 'i'
      IF (jformz==2) line(last2) = 'I'
      IF (jformz==3) line(last2) = ')'

      IF (line(kpt)==' ' .AND. line(kpt+1)=='+') THEN
        DO 60 j = kpt + 2, last2
          line(j-2) = line(j)
60      CONTINUE
        line(last2-1) = ' '
        line(last2) = ' '
        last2 = last2 - 2
      END IF

      ma(kptimu+2) = maim2
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmout
    SUBROUTINE zmpack(ma,mp)

!  MA is packed two base NDIG digits per word and returned in MP.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mp(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack
! ..
      CALL fmpack(ma,mp)
      CALL fmpack(ma(kptimu),mp(kptimp))
      RETURN
    END SUBROUTINE zmpack
    SUBROUTINE zmprnt(ma)

!  Print MA in base 10 format.

!  ZMPRNT can be called directly by the user for easy output
!  in M format.  MA is converted using ZMOUT and printed.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC int, log10, max, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz)
! ..
! .. Local Scalars ..
      INTEGER :: k, ksave, last1, last2, lb, lbz, nd, nexp
      CHARACTER (20) :: form
! ..
! .. External Subroutines ..
      EXTERNAL fmprnt, zmout
! ..
! .. Scalars in Common ..
      REAL (KIND(0.0D0)) :: dpmax
      REAL (KIND(0.0D0)) :: maxint, mbase, mexpov, mexpun, munkno, mxbase, mxexp, &
        mxexp2
      REAL :: runkno, spmax
      INTEGER :: intmax, iunkno, jform1, jform2, jformz, jprntz, kaccsw, &
        kdebug, keswch, kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, &
        ncall, ndg2mx, ndig, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mwa(lmwa)
      CHARACTER (1) :: cmbuff(lmbuff), cmbufz(lmbufz)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
      COMMON /zmbuff/cmbufz
      COMMON /zmuser/jformz, jprntz
! ..
      ksave = kflag
      nd = int(real(ndig)*log10(real(mbase))) + 1
      IF (nd<2) nd = 2
      nexp = int(2.0*log10(real(mxbase))) + 6
      lb = max(jform2+nexp,nd+nexp)

      IF (2*lb+7<=lmbufz .AND. jprntz==1) THEN
        lbz = 2*lb + 7
        CALL zmout(ma,cmbufz,lbz,last1,last2)
        WRITE (form,90000) kswide - 7
        WRITE (kw,form) (cmbufz(k),k=1,last2)
      ELSE
        CALL fmprnt(ma)
        CALL fmprnt(ma(kptimu))
      END IF
      kflag = ksave
      RETURN
90000 FORMAT (' (6X,',I3,'A1) ')
    END SUBROUTINE zmprnt
    SUBROUTINE zmpwr(ma,mb,mc)

!  MC = MA ** MB.

      IMPLICIT NONE

!             Scratch array usage during ZMPWR:   M01 - M06, MZ01 - MZ03

! .. Intrinsic Functions ..
      INTRINSIC abs, int, log, max, min, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz), mc(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccmb, maiz, marz, mbiz, mbrz, mtemp, mxsave
      REAL :: xval
      INTEGER :: iextra, intmb, j, jcos, jsin, jswap, k, kasave, kovun, &
        kradsv, kreslt, kwrnsv, ndsave
! ..
! .. External Functions ..
      LOGICAL, EXTERNAL :: fmcomp
! ..
! .. External Subroutines ..
      EXTERNAL fmabs, fmadd, fmcssn, fmdivi, fmexp, fmi2m, fmmi, fmmpy, &
        fmmpyd, fmmpyi, fmpi, fmrdc, zm2i2m, zmentr, zmeq, zmeq2, zmexit, &
        zmexp, zmi2m, zmipwr, zmln, zmmpy, zmmpyi, zmntr, zmsub, zmwarn
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
        mpisav(0:lunpck), mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      CALL zmentr('ZMPWR ',ma,mb,2,mc,kreslt,ndsave,mxsave,kasave,kovun)
      IF (kreslt/=0) RETURN
      marz = ma(0)
      mbrz = mb(0)
      maiz = ma(kptimu)
      mbiz = mb(kptimu)
      kaccsw = 0
      ndig = min(ndig+1,ndg2mx)

      CALL zmeq2(ma,ma,ndsave,ndig,1)
      CALL zmeq2(mb,mb,ndsave,ndig,1)

!             Check for special cases.

      IF (ma(2)==0 .AND. ma(kptimu+2)==0) THEN
        IF (mb(2)>0 .AND. mb(kptimu+2)==0) THEN
          CALL zmi2m(0,mz02)
          GO TO 60
        ELSE
          kflag = -4
          mz02(1) = munkno
          mz02(2) = 1
          mz02(kptimu+1) = munkno
          mz02(kptimu+2) = 1
          DO 10 j = 2, ndsave
            mz02(j+1) = 0
            mz02(kptimu+j+1) = 0
10        CONTINUE
          mz02(0) = nint(ndig*alogm2)
          mz02(kptimu) = nint(ndig*alogm2)
          GO TO 60
        END IF
      END IF
      IF (mb(kptimu+2)==0) THEN
        kwrnsv = kwarn
        kwarn = 0
        CALL fmmi(mb,intmb)
        kwarn = kwrnsv
        IF (kflag==0) THEN
          IF (ncall==1) THEN
            xval = abs(intmb) + 1
            k = int((1.5*log(xval))/alogmb+2.0)
            ndig = max(ndig+k,2)
            IF (ndig>ndg2mx) THEN
              kflag = -9
              CALL zmwarn
              mb(1) = munkno
              mb(2) = 1
              mb(kptimu+1) = munkno
              mb(kptimu+2) = 1
              DO 20 j = 2, ndsave
                mb(j+1) = 0
                mb(kptimu+j+1) = 0
20            CONTINUE
              ndig = ndsave
              IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
              ncall = ncall - 1
              RETURN
            END IF
            IF (mbase>=100*abs(ma(2)) .OR. mbase>=100*abs(ma(kptimu+2))) THEN
              ndig = min(ndig+1,ndg2mx)
            END IF
          END IF
          CALL zmeq2(ma,ma,ndsave,ndig,1)
          CALL zmipwr(ma,intmb,mz03)
          CALL zmeq(mz03,mz02)
          GO TO 60
        END IF
      END IF

!             Check for cases where ABS(MA) is very close to 1, and
!             avoid cancellation.

      CALL fmabs(ma,m03)
      CALL fmabs(ma(kptimu),m04)
      CALL fmi2m(1,m05)
      IF (fmcomp(m03,'EQ',m05) .AND. (m04(1)<=(-ndig) .OR. m04(2)==0)) THEN
        IF (ma(2)>0) THEN

!                 (1+c)**b = 1 + b*c + ...

          CALL zmi2m(1,mz02)
          CALL zmsub(ma,mz02,mz02)
          CALL zmmpy(mb,mz02,mz02)
          CALL fmadd(mz02,m05,mz02)
        ELSE

!                 (-1+c)**b = (-1)**b * (1 - b*c + ... )

          CALL zmi2m(-1,mz02)
          CALL zmsub(ma,mz02,mz02)
          CALL zmmpy(mb,mz02,mz02)
          CALL zmmpyi(mz02,-1,mz02)
          CALL fmadd(mz02,m05,mz02)
          kradsv = krad
          krad = 0
          IF (ma(kptimu+2)>=0) THEN
            CALL fmmpyi(mb,180,m06)
          ELSE
            CALL fmmpyi(mb,-180,m06)
          END IF
          CALL fmcssn(m06,mz03,mz03(kptimu))
          krad = kradsv
          CALL fmpi(m05)
          CALL fmmpy(m05,mb(kptimu),m05)
          IF (ma(kptimu+2)>=0) CALL fmmpyi(m05,-1,m05)
          CALL fmexp(m05,m05)
          CALL fmmpyd(m05,mz03,mz03(kptimu),mz03,mz03(kptimu))
          CALL zmmpy(mz02,mz03,mz02)
        END IF
        GO TO 60
      END IF
      IF (fmcomp(m04,'EQ',m05) .AND. (m03(1)<=(-ndig) .OR. m03(2)==0)) THEN
        IF (ma(kptimu+2)>0) THEN

!                 (i+c)**b = i**b * (1 - b*c*i - ... )

          CALL zm2i2m(0,1,mz02)
          CALL zmsub(ma,mz02,mz02)
          CALL zmmpy(mb,mz02,mz02)
          DO 30 j = 0, ndig + 1
            mtemp = mz02(j)
            mz02(j) = mz02(kptimu+j)
            mz02(kptimu+j) = mtemp
30        CONTINUE
          IF (mz02(kptimu+1)/=munkno) mz02(kptimu+2) = -mz02(kptimu+2)
          CALL fmadd(mz02,m05,mz02)
          kradsv = krad
          krad = 0
          CALL fmmpyi(mb,90,m06)
          CALL fmcssn(m06,mz03,mz03(kptimu))
          krad = kradsv
          CALL fmpi(m05)
          CALL fmmpy(m05,mb(kptimu),m05)
          CALL fmdivi(m05,-2,m05)
          CALL fmexp(m05,m05)
          CALL fmmpyd(m05,mz03,mz03(kptimu),mz03,mz03(kptimu))
          CALL zmmpy(mz02,mz03,mz02)
        ELSE

!                 (-i+c)**b = (-i)**b * (1 + b*c*i - ... )

          CALL zm2i2m(0,-1,mz02)
          CALL zmsub(ma,mz02,mz02)
          CALL zmmpy(mb,mz02,mz02)
          DO 40 j = 0, ndig + 1
            mtemp = mz02(j)
            mz02(j) = mz02(kptimu+j)
            mz02(kptimu+j) = mtemp
40        CONTINUE
          IF (mz02(1)/=munkno) mz02(2) = -mz02(2)
          CALL fmadd(mz02,m05,mz02)
          kradsv = krad
          krad = 0
          CALL fmmpyi(mb,-90,m06)
          CALL fmcssn(m06,mz03,mz03(kptimu))
          krad = kradsv
          CALL fmpi(m05)
          CALL fmmpy(m05,mb(kptimu),m05)
          CALL fmdivi(m05,2,m05)
          CALL fmexp(m05,m05)
          CALL fmmpyd(m05,mz03,mz03(kptimu),mz03,mz03(kptimu))
          CALL zmmpy(mz02,mz03,mz02)
        END IF
        GO TO 60
      END IF

      CALL zmln(ma,mz02)
      CALL zmmpy(mb,mz02,mz02)
      kwrnsv = kwarn
      kwarn = 0
      CALL fmrdc(mz02(kptimu),mz01,jsin,jcos,jswap)
      kwarn = kwrnsv
      IF (kflag==-9) THEN
        iextra = int(mz01(1))
      ELSE
        iextra = int(mz02(kptimu+1)-mz01(1))
      END IF
      IF (iextra>1) THEN
        ndig = ndig + iextra
        IF (ndig>ndg2mx) THEN
          kflag = -9
          CALL zmwarn
          mz02(1) = munkno
          mz02(2) = 1
          mz02(kptimu+1) = munkno
          mz02(kptimu+2) = 1
          DO 50 j = 2, ndsave
            mz02(j+1) = 0
            mz02(kptimu+j+1) = 0
50        CONTINUE
          ndig = ndig - iextra
          mz02(0) = nint(ndig*alogm2)
          mz02(kptimu) = nint(ndig*alogm2)
          GO TO 60
        END IF
        CALL zmeq2(ma,ma,ndsave,ndig,1)
        CALL zmeq2(mb,mb,ndsave,ndig,1)
        CALL zmln(ma,mz02)
        CALL zmmpy(mb,mz02,mz02)
      END IF

      CALL zmexp(mz02,mz02)

60    maccmb = mz02(0)
      ma(0) = marz
      mb(0) = mbrz
      mz02(0) = min(maccmb,marz,maiz,mbrz,mbiz)
      maccmb = mz02(kptimu)
      ma(kptimu) = maiz
      mb(kptimu) = mbiz
      mz02(kptimu) = min(maccmb,marz,maiz,mbrz,mbiz)
      CALL zmexit(mz02,mc,ndsave,mxsave,kasave,kovun,0)
      RETURN
    END SUBROUTINE zmpwr
    SUBROUTINE zmread(kread,ma)

!  Read MA on unit KREAD.  Multi-line numbers will have '&' as the
!  last nonblank character on all but the last line.  Only one
!  number is allowed on the line(s).

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kread
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz)
! ..
! .. Local Scalars ..
      INTEGER :: j, lb
! ..
! .. Local Arrays ..
      CHARACTER (1) :: line(80)
! ..
! .. External Subroutines ..
      EXTERNAL zminp, zmwarn
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
      CHARACTER (1) :: cmbuff(lmbuff), cmbufz(lmbufz)
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
      COMMON /zmbuff/cmbufz
! ..
      ncall = ncall + 1
      namest(ncall) = 'ZMREAD'
      lb = 0

10    READ (kread,90000,err=30,end=30) line

!             Scan the line and look for '&'

      DO 20 j = 1, 80
        IF (line(j)=='&') GO TO 10
        IF (line(j)/=' ') THEN
          lb = lb + 1
          IF (lb>lmbufz) THEN
            kflag = -8
            GO TO 40
          END IF
          cmbufz(lb) = line(j)
        END IF
20    CONTINUE

      CALL zminp(cmbufz,ma,1,lb)

      ncall = ncall - 1
      RETURN

!             If there is an error, return UNKNOWN.

30    kflag = -4
40    CALL zmwarn
      ma(1) = munkno
      ma(2) = 1
      ma(kptimu+1) = munkno
      ma(kptimu+2) = 1
      ma(0) = nint(ndig*alogm2)
      ma(kptimu) = nint(ndig*alogm2)
      DO 50 j = 2, ndig
        ma(j+1) = 0
        ma(kptimu+j+1) = 0
50    CONTINUE
      ncall = ncall - 1
      RETURN
90000 FORMAT (80A1)
    END SUBROUTINE zmread
    SUBROUTINE zmreal(ma,mbfm)

!  MBFM = REAL(MA)

!  MA is a complex ZM number, MBFM is a real FM number.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mbfm(0:lunpck)
! ..
! .. External Subroutines ..
      EXTERNAL fmeq, fmntr, zmntr
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
      kflag = 0
      ncall = ncall + 1
      namest(ncall) = 'ZMREAL'
      IF (ntrace/=0) CALL zmntr(2,ma,ma,1)

      CALL fmeq(ma,mbfm)

      IF (ntrace/=0) CALL fmntr(1,mbfm,mbfm,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmreal
    SUBROUTINE zmrpwr(ma,ival,jval,mb)

!  MB = MA ** (IVAL/JVAL)

!  Raise a ZM number to a rational power.

      IMPLICIT NONE

!             Scratch array usage during ZMRPWR:  M01 - M03, MZ01 - MZ04

! .. Intrinsic Functions ..
      INTRINSIC abs, atan2, cos, dble, int, log, max, min, nint, real, sin
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival, jval
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: ar, br, f, theta, x
      REAL (KIND(0.0D0)) :: ma2, maccmb, maiz, marz, mr1, mxsave
      REAL :: xval
      INTEGER :: ijsign, invert, ival2, j, jval2, k, kasave, kovun, kst, l, &
        lval, ndsave
! ..
! .. Local Arrays ..
      INTEGER :: nstack(19)
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmcons, fmdig, fmdiv, fmdpm, fmgcdi, fmm2dp, fmmpy, &
        fmntri, fmsqr, fmsqrt, zmadd, zmdiv, zmdivi, zmeq2, zmexit, zmi2m, &
        zmipwr, zmmpyi, zmntr, zmsqrt, zmwarn
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
        mpisav(0:lunpck), mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      ncall = ncall + 1
      namest(ncall) = 'ZMRPWR'
      ndsave = ndig
      IF (ntrace/=0) THEN
        CALL zmntr(2,ma,ma,1)
        CALL fmntri(2,ival,0)
        CALL fmntri(2,jval,0)
      END IF
      marz = ma(0)
      maiz = ma(kptimu)
      kovun = 0
      IF (ma(1)==mexpov .OR. ma(1)==mexpun .OR. ma(kptimu+1)==mexpov .OR. &
        ma(kptimu+1)==mexpun) kovun = 1

      IF (mblogs/=mbase) CALL fmcons
      kflag = 0
      ijsign = 1
      ival2 = abs(ival)
      jval2 = abs(jval)
      IF (ival>0 .AND. jval<0) ijsign = -1
      IF (ival<0 .AND. jval>0) ijsign = -1
      IF (ival2>0 .AND. jval2>0) CALL fmgcdi(ival2,jval2)

!             Check for special cases.

      IF (ma(1)==munkno .OR. ma(kptimu+1)==munkno .OR. (ijsign<=0 .AND. ma( &
          2)==0 .AND. ma(kptimu+2)==0) .OR. jval==0) THEN
        ma2 = ma(2)
        mb(0) = nint(ndig*alogm2)
        mb(1) = munkno
        mb(2) = 1
        mb(kptimu) = nint(ndig*alogm2)
        mb(kptimu+1) = munkno
        mb(kptimu+2) = 1
        DO 10 j = 2, ndsave
          mb(j+1) = 0
          mb(kptimu+j+1) = 0
10      CONTINUE
        kflag = -4
        IF (ival<=0 .AND. ma2==0) CALL zmwarn
        IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
        ncall = ncall - 1
        RETURN
      END IF

      IF (ival==0) THEN
        CALL zmi2m(1,mb)
        IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
        ncall = ncall - 1
        RETURN
      END IF

!             Increase the working precision.

      IF (ncall==1) THEN
        xval = max(abs(ival),abs(jval)) + 1
        k = int((5.0*real(dlogtn)+log(xval))/alogmb+2.0)
        ndig = max(ndig+k,2)
      ELSE
        xval = max(abs(ival),abs(jval)) + 1
        k = int(log(xval)/alogmb+1.0)
        ndig = ndig + k
      END IF
      IF (ndig>ndg2mx) THEN
        kflag = -9
        CALL zmwarn
        mb(1) = munkno
        mb(2) = 1
        mb(kptimu+1) = munkno
        mb(kptimu+2) = 1
        DO 20 j = 2, ndsave
          mb(j+1) = 0
          mb(kptimu+j+1) = 0
20      CONTINUE
        ndig = ndsave
        mb(0) = nint(ndig*alogm2)
        mb(kptimu) = nint(ndig*alogm2)
        ndig = ndsave
        IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
        ncall = ncall - 1
        RETURN
      END IF
      IF (mbase>=100*abs(ma(2)) .OR. mbase>=100*abs(ma(kptimu+2))) THEN
        ndig = min(ndig+1,ndg2mx)
      END IF
      kasave = kaccsw
      mxsave = mxexp
      mxexp = mxexp2

      CALL zmeq2(ma,mz04,ndsave,ndig,0)
      IF (ival2==1 .AND. jval2==2) THEN
        CALL zmsqrt(mz04,mb)
        GO TO 50
      END IF

!             Generate the first approximation to MA**(1/JVAL2).

      CALL zmi2m(0,mb)
      CALL fmdig(nstack,kst)
      ndig = nstack(1)
      CALL fmsqr(mz04,mz03)
      CALL fmsqr(mz04(kptimu),m03)
      CALL fmadd(mz03,m03,mz03)
      CALL fmsqrt(mz03,mz03)
      IF (mz03(1)>=mexpov) THEN
        kflag = -4
        CALL zmwarn
        mb(1) = munkno
        mb(2) = 1
        mb(kptimu+1) = munkno
        mb(kptimu+2) = 1
        DO 30 j = 2, ndsave
          mb(j+1) = 0
          mb(kptimu+j+1) = 0
30      CONTINUE
        ndig = ndsave
        mb(0) = nint(ndig*alogm2)
        mb(kptimu) = nint(ndig*alogm2)
        mxexp = mxsave
        kaccsw = kasave
        ndig = ndsave
        IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
        ncall = ncall - 1
        RETURN
      END IF

!             Invert MA if ABS(MA) > 1 and IVAL or JVAL is large.

      invert = 0
      IF (ival>5 .OR. jval>5) THEN
        IF (mz03(1)>0) THEN
          invert = 1
          ndig = nstack(kst)
          CALL zmi2m(1,mb)
          CALL zmdiv(mb,mz04,mz04)
          ndig = nstack(1)
          CALL fmdiv(mb,mz03,mz03)
        END IF
      END IF

      CALL fmdiv(mz04,mz03,m03)
      CALL fmm2dp(m03,ar)
      CALL fmdiv(mz04(kptimu),mz03,m03)
      CALL fmm2dp(m03,br)
      mr1 = mz03(1)
      mz03(1) = 0
      CALL fmm2dp(mz03,x)
      l = int(mr1/jval2)
      f = mr1/dble(jval2) - l
      x = x**(1.0D0/jval2)*dble(mbase)**f
      CALL fmdpm(x,m03)
      m03(1) = m03(1) + l

      theta = atan2(br,ar)
      x = cos(theta/jval2)
      CALL fmdpm(x,mb)
      x = sin(theta/jval2)
      CALL fmdpm(x,mb(kptimu))
      CALL fmmpy(m03,mb,mb)
      CALL fmmpy(m03,mb(kptimu),mb(kptimu))

!             Newton iteration.

      DO 40 j = 1, kst
        ndig = nstack(j)
        IF (j<kst) ndig = ndig + 1
        lval = jval2 - 1
        CALL zmipwr(mb,lval,mz03)
        CALL zmdiv(mz04,mz03,mz03)
        CALL zmmpyi(mb,lval,mb)
        CALL zmadd(mb,mz03,mb)
        CALL zmdivi(mb,jval2,mb)
40    CONTINUE

      CALL zmipwr(mb,ijsign*ival2,mb)
      IF (invert==1) THEN
        CALL zmi2m(1,mz03)
        CALL zmdiv(mz03,mb,mb)
      END IF

!             Round the result and return.

50    maccmb = mb(0)
      ma(0) = marz
      mb(0) = min(maccmb,marz,maiz)
      maccmb = mb(kptimu)
      ma(kptimu) = maiz
      mb(kptimu) = min(maccmb,marz,maiz)
      CALL zmexit(mb,mb,ndsave,mxsave,kasave,kovun,1)
      RETURN
    END SUBROUTINE zmrpwr
    SUBROUTINE zmrslt(mc,kreslt)

!  Handle results that are special cases, such as overflow,
!  underflow, and unknown.

!  MC is the result that is returned

!  KRESLT is the result code.  Result codes handled here:

!   0 - Perform the normal operation
!  12 - The result is 'UNKNOWN'

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kreslt
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: mc(0:lunpkz)
! ..
! .. Local Scalars ..
      INTEGER :: kfsave
! ..
! .. External Subroutines ..
      EXTERNAL zmi2m
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

      IF (kreslt==12 .OR. kreslt<0 .OR. kreslt>15) THEN
        CALL zmi2m(0,mc)
        mc(1) = munkno
        mc(2) = 1
        mc(0) = nint(ndig*alogm2)
        mc(kptimu+1) = munkno
        mc(kptimu+2) = 1
        mc(kptimu) = nint(ndig*alogm2)
        kflag = kfsave
        RETURN
      END IF

      RETURN
    END SUBROUTINE zmrslt
    SUBROUTINE zmsin(ma,mb)

!  MB = SIN(MA).

      IMPLICIT NONE

!             Scratch array usage during ZMSIN:   M01 - M06, MZ01 - MZ03

! .. Intrinsic Functions ..
      INTRINSIC min
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccmb, maiz, marz, mxsave
      INTEGER :: kasave, kovun, kreslt, krsave, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmchsh, fmcssn, fmi2m, fmmpy, fmsin, fmsinh, zmentr, zmeq, &
        zmeq2, zmexit, zmi2m
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
        mpisav(0:lunpck), mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      CALL zmentr('ZMSIN ',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
      IF (kreslt/=0) RETURN
      marz = ma(0)
      maiz = ma(kptimu)
      kaccsw = 0
      krsave = krad
      krad = 1

      CALL zmeq2(ma,ma,ndsave,ndig,1)

!             Check for special cases.

      IF (ma(2)==0 .AND. ma(kptimu+2)==0) THEN
        CALL zmi2m(0,mz01)
        GO TO 10
      ELSE IF (ma(1)<(-ndig) .AND. ma(kptimu+1)<(-ndig)) THEN
        CALL zmeq(ma,mz01)
        GO TO 10
      ELSE IF (ma(kptimu+2)==0) THEN
        CALL fmsin(ma,mz01)
        CALL fmi2m(0,mz01(kptimu))
        GO TO 10
      ELSE IF (ma(2)==0) THEN
        CALL fmsinh(ma(kptimu),mz01(kptimu))
        CALL fmi2m(0,mz01)
        GO TO 10
      END IF

!             Find COS(REAL(MA)) and SIN(REAL(MA)).

      CALL fmcssn(ma,mz01(kptimu),mz01)

!             Find COSH(IMAG(MA)) and SINH(IMAG(MA)).

      CALL fmchsh(ma(kptimu),m05,m06)

!             SIN(MA) =  SIN(REAL(MA))*COSH(IMAG(MA)) +
!                        COS(REAL(MA))*SINH(IMAG(MA)) i

      CALL fmmpy(mz01,m05,mz01)
      CALL fmmpy(mz01(kptimu),m06,mz01(kptimu))

10    maccmb = mz01(0)
      ma(0) = marz
      mz01(0) = min(maccmb,marz,maiz)
      maccmb = mz01(kptimu)
      ma(kptimu) = maiz
      mz01(kptimu) = min(maccmb,marz,maiz)
      CALL zmexit(mz01,mb,ndsave,mxsave,kasave,kovun,0)
      krad = krsave
      RETURN
    END SUBROUTINE zmsin
    SUBROUTINE zmsinh(ma,mb)

!  MB = SINH(MA).

      IMPLICIT NONE

!             Scratch array usage during ZMSINH:  M01 - M06, MZ01 - MZ03

! .. Intrinsic Functions ..
      INTRINSIC min
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccmb, maiz, marz, mxsave
      INTEGER :: kasave, kovun, kreslt, krsave, ndsave
! ..
! .. External Subroutines ..
      EXTERNAL fmchsh, fmcssn, fmi2m, fmmpy, fmsin, fmsinh, zmentr, zmeq, &
        zmeq2, zmexit, zmi2m
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
        mpisav(0:lunpck), mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      CALL zmentr('ZMSINH',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
      IF (kreslt/=0) RETURN
      marz = ma(0)
      maiz = ma(kptimu)
      kaccsw = 0
      krsave = krad
      krad = 1

      CALL zmeq2(ma,ma,ndsave,ndig,1)

!             Check for special cases.

      IF (ma(2)==0 .AND. ma(kptimu+2)==0) THEN
        CALL zmi2m(0,mz01)
        GO TO 10
      ELSE IF (ma(1)<(-ndig) .AND. ma(kptimu+1)<(-ndig)) THEN
        CALL zmeq(ma,mz01)
        GO TO 10
      ELSE IF (ma(2)==0) THEN
        CALL fmsin(ma(kptimu),mz01(kptimu))
        CALL fmi2m(0,mz01)
        GO TO 10
      ELSE IF (ma(kptimu+2)==0) THEN
        CALL fmsinh(ma,mz01)
        CALL fmi2m(0,mz01(kptimu))
        GO TO 10
      END IF

!             Find SIN(IMAG(MA)) and COS(IMAG(MA)).

      CALL fmcssn(ma(kptimu),mz01,mz01(kptimu))

!             Find SINH(REAL(MA)) and COSH(REAL(MA)).

      CALL fmchsh(ma,m05,m06)

!             SINH(MA) =  SINH(REAL(MA))*COS(IMAG(MA)) +
!                         COSH(REAL(MA))*SIN(IMAG(MA)) i

      CALL fmmpy(mz01,m06,mz01)
      CALL fmmpy(mz01(kptimu),m05,mz01(kptimu))

10    maccmb = mz01(0)
      ma(0) = marz
      mz01(0) = min(maccmb,marz,maiz)
      maccmb = mz01(kptimu)
      ma(kptimu) = maiz
      mz01(kptimu) = min(maccmb,marz,maiz)
      CALL zmexit(mz01,mb,ndsave,mxsave,kasave,kovun,0)
      krad = krsave
      RETURN
    END SUBROUTINE zmsinh
    SUBROUTINE zmsqr(ma,mb)

!  MB = MA * MA

      IMPLICIT NONE

!             Scratch array usage during ZMSQR:   M01 - M03, MZ01

! .. Intrinsic Functions ..
      INTRINSIC abs, max, min, nint
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccmb, maiz, marz, mxsave, mzero
      INTEGER :: j, kasave, kovun, kreslt, kwrnsv, ndgsv2, ndsave, ntrsav
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmi2m, fmmpy, fmsqr, fmsub, zmentr, zmeq2, zmntr, &
        zmrslt, zmwarn
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
        mpisav(0:lunpck), mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      IF (abs(ma(1))>mexpab .OR. abs(ma(kptimu+1))>mexpab .OR. kdebug>=1) THEN
        CALL zmentr('ZMSQR ',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        IF (ntrace/=0) THEN
          namest(ncall) = 'ZMSQR '
          CALL zmntr(2,ma,ma,1)
        END IF
        ndsave = ndig
        IF (ncall==1) THEN
          ndig = max(ndig+ngrd52,2)
          IF (ndig>ndg2mx) THEN
            namest(ncall) = 'ZMSQR '
            kflag = -9
            CALL zmwarn
            kreslt = 12
            ndig = ndsave
            CALL zmrslt(mb,kreslt)
            IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
            ncall = ncall - 1
            RETURN
          END IF
          IF (mbase>=100*abs(ma(2)) .OR. mbase>=100*abs(ma(kptimu+2))) THEN
            ndig = min(ndig+1,ndg2mx)
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 0
        mxsave = mxexp
        mxexp = mxexp2
        kovun = 0
      END IF

      marz = ma(0)
      maiz = ma(kptimu)
      mzero = 0
      ntrsav = ntrace
      ntrace = 0
      kwrnsv = kwarn
      kwarn = 0

      DO 10 j = ndsave + 2, ndig + 1
        ma(j) = mzero
        ma(kptimu+j) = mzero
10    CONTINUE
      IF (ncall==1) THEN
        ma(0) = nint(ndig*alogm2)
        ma(kptimu) = ma(0)
      END IF

!             Check for special cases.

      IF (ma(kptimu+2)==0) THEN
        CALL fmsqr(ma,mz01)
        CALL fmi2m(0,mz01(kptimu))
      ELSE IF (ma(2)==0) THEN
        CALL fmsqr(ma(kptimu),mz01)
        IF (mz01(1)/=munkno) mz01(2) = -mz01(2)
        CALL fmi2m(0,mz01(kptimu))
      ELSE
        CALL fmadd(ma,ma(kptimu),m02)
        CALL fmsub(ma,ma(kptimu),m03)
        CALL fmmpy(m02,m03,mz01)
        CALL fmmpy(ma,ma(kptimu),m03)
        CALL fmadd(m03,m03,mz01(kptimu))
      END IF

      mxexp = mxsave
      ntrace = ntrsav
      ndgsv2 = ndig
      ndig = ndsave
      kwarn = kwrnsv
      IF (ncall==1) THEN
        ma(0) = marz
        ma(kptimu) = maiz
      END IF
      maccmb = mz01(0)
      ma(0) = marz
      mz01(0) = min(maccmb,marz,maiz)
      maccmb = mz01(kptimu)
      ma(kptimu) = maiz
      mz01(kptimu) = min(maccmb,marz,maiz)
      kaccsw = kasave
      CALL zmeq2(mz01,mb,ndgsv2,ndsave,0)
      IF (mb(1)>=mexpov .OR. mb(1)<=-mexpov .OR. mb(kptimu+1)>=mexpov .OR. &
          mb(kptimu+1)<=-mexpov) THEN
        IF (mb(1)==munkno .OR. mb(kptimu+1)==munkno) THEN
          kflag = -4
        ELSE IF (mb(1)==mexpov .OR. mb(kptimu+1)==mexpov) THEN
          kflag = -5
        ELSE IF (mb(1)==mexpun .OR. mb(kptimu+1)==mexpun) THEN
          kflag = -6
        END IF
        IF ((mb(1)==munkno) .OR. (mb(kptimu+1)==munkno) .OR. (mb(1)==mexpun &
            .AND. kovun==0) .OR. (mb(kptimu+1)==mexpun .AND. kovun==0) .OR. ( &
            mb(1)==mexpov .AND. kovun==0) .OR. (mb(kptimu+ &
            1)==mexpov .AND. kovun==0)) THEN
          namest(ncall) = 'ZMSQR '
          CALL zmwarn
        END IF
      END IF
      IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmsqr
    SUBROUTINE zmsqrt(ma,mb)

!  MB = SQRT(MA).  Principal Square Root.

      IMPLICIT NONE

!             Scratch array usage during ZMSQRT:   M01 - M03, MZ01

! .. Intrinsic Functions ..
      INTRINSIC abs, max, min
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccmb, maiz, marz, mxsave
      INTEGER :: kasave, kovun, kreslt, kwrnsv, ndsave, ntrsav
! ..
! .. External Subroutines ..
      EXTERNAL fmabs, fmadd, fmdiv, fmdivi, fmeq, fmsqr, fmsqrt, zmentr, &
        zmeq2, zmi2m, zmntr, zmrslt, zmwarn
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
        mpisav(0:lunpck), mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      IF (abs(ma(1))>mexpab .OR. abs(ma(kptimu+1))>mexpab .OR. kdebug>=1) THEN
        CALL zmentr('ZMSQRT',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
      ELSE
        ncall = ncall + 1
        IF (ntrace/=0) THEN
          namest(ncall) = 'ZMSQRT'
          CALL zmntr(2,ma,ma,1)
        END IF
        ndsave = ndig
        IF (ncall==1) THEN
          ndig = max(ndig+ngrd52,2)
          IF (ndig>ndg2mx) THEN
            namest(ncall) = 'ZMSQRT'
            kflag = -9
            CALL zmwarn
            kreslt = 12
            ndig = ndsave
            CALL zmrslt(mb,kreslt)
            IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
            ncall = ncall - 1
            RETURN
          END IF
          IF (mbase>=100*abs(ma(2)) .OR. mbase>=100*abs(ma(kptimu+2))) THEN
            ndig = min(ndig+1,ndg2mx)
          END IF
        END IF
        kasave = kaccsw
        kaccsw = 0
        mxsave = mxexp
        mxexp = mxexp2
        kovun = 0
      END IF

      ntrsav = ntrace
      ntrace = 0
      kwrnsv = kwarn
      kwarn = 0
      marz = ma(0)
      maiz = ma(kptimu)

      CALL zmeq2(ma,ma,ndsave,ndig,1)

!             Check for special cases.

      IF (ma(2)==0 .AND. ma(kptimu+2)==0) THEN
        CALL zmi2m(0,mz01)
        GO TO 10
      ELSE IF (ma(2)==0) THEN
        CALL fmabs(ma(kptimu),m01)
        CALL fmdivi(m01,2,m03)
        CALL fmsqrt(m03,m03)
      ELSE IF (ma(kptimu+2)==0) THEN
        CALL fmabs(ma,m03)
        CALL fmsqrt(m03,m03)
      ELSE
        CALL fmsqr(ma,m01)
        CALL fmsqr(ma(kptimu),m02)
        CALL fmadd(m01,m02,m03)
        CALL fmsqrt(m03,m03)
        CALL fmabs(ma,m02)
        CALL fmadd(m02,m03,m03)
        CALL fmdivi(m03,2,m03)
        CALL fmsqrt(m03,m03)
      END IF

      CALL fmadd(m03,m03,m02)
      IF (ma(2)>=0) THEN
        CALL fmdiv(ma(kptimu),m02,mz01(kptimu))
        CALL fmeq(m03,mz01)
      ELSE
        IF (ma(kptimu+2)>=0) THEN
          CALL fmdiv(ma(kptimu),m02,mz01)
          CALL fmeq(m03,mz01(kptimu))
        ELSE
          CALL fmdiv(ma(kptimu),m02,mz01)
          CALL fmeq(m03,mz01(kptimu))
          IF (mz01(1)/=munkno) mz01(2) = -mz01(2)
          IF (mz01(kptimu+1)/=munkno) mz01(kptimu+2) = -mz01(kptimu+2)
        END IF
      END IF

10    mxexp = mxsave
      maccmb = mz01(0)
      ma(0) = marz
      mz01(0) = min(maccmb,marz,maiz)
      maccmb = mz01(kptimu)
      ma(kptimu) = maiz
      mz01(kptimu) = min(maccmb,marz,maiz)
      kaccsw = kasave
      CALL zmeq2(mz01,mb,ndig,ndsave,0)

      IF (mb(1)==munkno .OR. mb(kptimu+1)==munkno) THEN
        kflag = -4
      ELSE IF (mb(1)==mexpov .OR. mb(kptimu+1)==mexpov) THEN
        kflag = -5
      ELSE IF (mb(1)==mexpun .OR. mb(kptimu+1)==mexpun) THEN
        kflag = -6
      END IF
      ntrace = ntrsav
      ndig = ndsave
      kwarn = kwrnsv
      IF ((mb(1)==munkno) .OR. (mb(kptimu+1)==munkno) .OR. (mb(1)==mexpun &
          .AND. kovun==0) .OR. (mb(kptimu+1)==mexpun .AND. kovun==0) .OR. (mb( &
          1)==mexpov .AND. kovun==0) .OR. (mb(kptimu+ &
          1)==mexpov .AND. kovun==0)) THEN
        namest(ncall) = 'ZMSQRT'
        CALL zmwarn
      END IF
      IF (ntrace/=0) CALL zmntr(1,mb,mb,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmsqrt
    SUBROUTINE zmst2m(string,ma)

!  MA = STRING

!  Convert a character string to FM format.
!  This is often more convenient than using ZMINP, which converts an
!  array of CHARACTER*1 values.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC len
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (*) :: string
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz)
! ..
! .. Local Scalars ..
      INTEGER :: j, lb
! ..
! .. External Subroutines ..
      EXTERNAL zminp
! ..
! .. Scalars in Common ..
      REAL :: alogm2, alogmb, alogmt, alogmx, runkno, spmax
      REAL (KIND(0.0D0)) :: dlogeb, dlogmb, dlogpi, dlogtn, dlogtp, &
        dlogtw, dpeps, dpmax, dppi
      REAL (KIND(0.0D0)) :: maxint, mbase, mblogs, mbse, mbslb, mbsli,  &
        mbspi, mexpab, mexpov, mexpun, munkno, mxbase, mxexp, mxexp2
      INTEGER :: intmax, iunkno, jform1, jform2, jformz, jprntz, kaccsw, &
        kdebug, keswch, kflag, krad, kround, ksub, kswide, kw, kwarn, lvltrc, &
        ncall, ndg2mx, ndig, ndige, ndiglb, ndigli, ndigpi, ngrd21, ngrd22, &
        ngrd52, ntrace
      CHARACTER (1) :: cmchar
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mesav(0:lunpck), mlbsav(0:lunpck), mln1(0:lunpck), &
        mln2(0:lunpck), mln3(0:lunpck), mln4(0:lunpck), mpisav(0:lunpck), &
        mwa(lmwa)
      INTEGER :: khasht(lhash1:lhash2), khashv(lhash1:lhash2)
      CHARACTER (1) :: cmbuff(lmbuff), cmbufz(lmbufz)
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
      COMMON /zmbuff/cmbufz
      COMMON /zmuser/jformz, jprntz
! ..
      ncall = ncall + 1
      namest(ncall) = 'ZMST2M'
      lb = len(string)

      DO 10 j = 1, lb
        cmbufz(j) = string(j:j)
10    CONTINUE

      CALL zminp(cmbufz,ma,1,lb)

      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmst2m
    SUBROUTINE zmsub(ma,mb,mc)

!  MC = MA - MB

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC abs
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz), mc(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: mxsave
      INTEGER :: kasave, kf1, kovun, kreslt, kwrnsv, ndsave, ntrsav
! ..
! .. External Subroutines ..
      EXTERNAL fmsub, zmentr, zmntr, zmwarn
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
      IF (abs(ma(1))>mexpab .OR. abs(ma(kptimu+1))>mexpab .OR. abs(mb(1))> &
          mexpab .OR. abs(mb(kptimu+1))>mexpab .OR. kdebug>=1) THEN
        CALL zmentr('ZMSUB ',ma,mb,2,mc,kreslt,ndsave,mxsave,kasave,kovun)
        IF (kreslt/=0) RETURN
        ndig = ndsave
        mxexp = mxsave
        kaccsw = kasave
      ELSE
        ncall = ncall + 1
        IF (ntrace/=0) THEN
          namest(ncall) = 'ZMSUB '
          CALL zmntr(2,ma,mb,2)
        END IF
        kovun = 0
      END IF

!             Force FMSUB to use more guard digits for user calls.

      ncall = ncall - 1
      ntrsav = ntrace
      ntrace = 0
      kwrnsv = kwarn
      kwarn = 0

      CALL fmsub(ma,mb,mc)
      kf1 = kflag
      CALL fmsub(ma(kptimu),mb(kptimu),mc(kptimu))

      ntrace = ntrsav
      kwarn = kwrnsv
      ncall = ncall + 1
      IF (ntrace/=0) namest(ncall) = 'ZMSUB '
      IF (kflag==1) kflag = kf1

      IF (mc(1)==munkno .OR. mc(kptimu+1)==munkno) THEN
        kflag = -4
      ELSE IF (mc(1)==mexpov .OR. mc(kptimu+1)==mexpov) THEN
        kflag = -5
      ELSE IF (mc(1)==mexpun .OR. mc(kptimu+1)==mexpun) THEN
        kflag = -6
      END IF
      IF ((mc(1)==munkno) .OR. (mc(kptimu+1)==munkno) .OR. (mc(1)==mexpun &
          .AND. kovun==0) .OR. (mc(kptimu+1)==mexpun .AND. kovun==0) .OR. (mc( &
          1)==mexpov .AND. kovun==0) .OR. (mc(kptimu+ &
          1)==mexpov .AND. kovun==0)) THEN
        namest(ncall) = 'ZMSUB '
        CALL zmwarn
      END IF
      IF (ntrace/=0) CALL zmntr(1,mc,mc,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmsub
    SUBROUTINE zmtan(ma,mb)

!  MB = TAN(MA).

      IMPLICIT NONE

!             Scratch array usage during ZMTAN:   M01 - M06, MZ01 - MZ03

! .. Intrinsic Functions ..
      INTRINSIC int, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccmb, maiz, marz, mxsave
      INTEGER :: iextra, j, kasave, kovun, kreslt, krsave, ndsave, ngoal
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmchsh, fmcssn, fmdiv, fmdivd, fmi2m, fmim, fmtan, &
        fmtanh, zmentr, zmeq, zmeq2, zmexit, zmi2m, zmwarn
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
        mpisav(0:lunpck), mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      CALL zmentr('ZMTAN ',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
      IF (kreslt/=0) RETURN
      marz = ma(0)
      maiz = ma(kptimu)
      krsave = krad
      krad = 1

10    CALL zmeq2(ma,ma,ndsave,ndig,1)

!             Check for special cases.

      IF (ma(2)==0 .AND. ma(kptimu+2)==0) THEN
        CALL zmi2m(0,mz01)
        GO TO 20
      ELSE IF (ma(1)<(-ndig) .AND. ma(kptimu+1)<(-ndig)) THEN
        CALL zmeq(ma,mz01)
        GO TO 20
      ELSE IF (ma(kptimu+2)==0) THEN
        CALL fmtan(ma,mz01)
        CALL fmi2m(0,mz01(kptimu))
        GO TO 20
      ELSE IF (ma(2)==0) THEN
        CALL fmtanh(ma(kptimu),mz01(kptimu))
        CALL fmi2m(0,mz01)
        GO TO 20
      END IF

!             Find SIN(2*REAL(MA)) and COS(2*REAL(MA)).

      CALL fmadd(ma,ma,mz01)
      CALL fmcssn(mz01,mz01(kptimu),mz01)

!             Find SINH(2*IMAG(MA)) and COSH(2*IMAG(MA)).

      CALL fmadd(ma(kptimu),ma(kptimu),m06)
      CALL fmchsh(m06,m05,m06)

!             TAN(MA) =  SIN(2*REAL(MA)) /
!                        (COS(2*REAL(MA))+COSH(2*IMAG(MA)) +
!                        SINH(2*IMAG(MA)) /
!                        (COS(2*REAL(MA))+COSH(2*IMAG(MA)) i

      CALL fmadd(mz01(kptimu),m05,m05)
      IF (m05(2)==0) THEN
        mz01(0) = 0
        ngoal = int(real(ndsave)*alogm2) + 7
        GO TO 30
      ELSE IF (m05(1)==mexpov) THEN
        CALL fmdiv(mz01,m05,mz01)
        CALL fmim(1,mz01(kptimu))
        IF (m06(2)<0) mz01(kptimu+2) = -mz01(kptimu+2)
      ELSE
        CALL fmdivd(mz01,m06,m05,mz01,mz01(kptimu))
      END IF

!             Check for too much cancellation.

20    IF (ncall<=1) THEN
        ngoal = int(real(ndsave)*alogm2) + 7
      ELSE
        ngoal = int(-mxexp2)
      END IF
30    IF (mz01(0)<=ngoal .OR. mz01(kptimu)<=ngoal) THEN
        iextra = int(real(max(ngoal-mz01(0),ngoal-mz01(kptimu)))/alogm2+23.03/ &
          alogmb) + 1
        ndig = ndig + iextra
        IF (ndig>ndg2mx) THEN
          kflag = -9
          CALL zmwarn
          mz01(1) = munkno
          mz01(2) = 1
          mz01(kptimu+1) = munkno
          mz01(kptimu+2) = 1
          DO 40 j = 2, ndsave
            mz01(j+1) = 0
            mz01(kptimu+j+1) = 0
40        CONTINUE
          ndig = ndig - iextra
          mz01(0) = nint(ndig*alogm2)
          mz01(kptimu) = nint(ndig*alogm2)
          GO TO 50
        END IF
        GO TO 10
      END IF

50    maccmb = mz01(0)
      ma(0) = marz
      mz01(0) = min(maccmb,marz,maiz)
      maccmb = mz01(kptimu)
      ma(kptimu) = maiz
      mz01(kptimu) = min(maccmb,marz,maiz)
      CALL zmexit(mz01,mb,ndsave,mxsave,kasave,kovun,0)
      krad = krsave
      RETURN
    END SUBROUTINE zmtan
    SUBROUTINE zmtanh(ma,mb)

!  MB = TANH(MA).

      IMPLICIT NONE

!             Scratch array usage during ZMTANH:  M01 - M06, MZ01 - MZ03

! .. Intrinsic Functions ..
      INTRINSIC int, max, min, nint, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: lhash1 = 0, lhash2 = 256, nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mb(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: maccmb, maiz, marz, mxsave
      INTEGER :: iextra, j, kasave, kovun, kreslt, krsave, ndsave, ngoal
! ..
! .. External Subroutines ..
      EXTERNAL fmadd, fmchsh, fmcssn, fmdiv, fmdivd, fmi2m, fmim, fmtan, &
        fmtanh, zmentr, zmeq, zmeq2, zmexit, zmi2m, zmwarn
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
        mpisav(0:lunpck), mwa(lmwa), mz01(0:lunpkz), mz02(0:lunpkz), &
        mz03(0:lunpkz), mz04(0:lunpkz)
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
      COMMON /zm1/mz01, mz02, mz03, mz04
! ..
      CALL zmentr('ZMTANH',ma,ma,1,mb,kreslt,ndsave,mxsave,kasave,kovun)
      IF (kreslt/=0) RETURN
      marz = ma(0)
      maiz = ma(kptimu)
      krsave = krad
      krad = 1

10    CALL zmeq2(ma,ma,ndsave,ndig,1)

!             Check for special cases.

      IF (ma(2)==0 .AND. ma(kptimu+2)==0) THEN
        CALL zmi2m(0,mz01)
        GO TO 20
      ELSE IF (ma(1)<(-ndig) .AND. ma(kptimu+1)<(-ndig)) THEN
        CALL zmeq(ma,mz01)
        GO TO 20
      ELSE IF (ma(2)==0) THEN
        CALL fmtan(ma(kptimu),mz01(kptimu))
        CALL fmi2m(0,mz01)
        GO TO 20
      ELSE IF (ma(kptimu+2)==0) THEN
        CALL fmtanh(ma,mz01)
        CALL fmi2m(0,mz01(kptimu))
        GO TO 20
      END IF

!             Find SIN(2*IMAG(MA)) and COS(2*IMAG(MA)).

      CALL fmadd(ma(kptimu),ma(kptimu),mz01)
      CALL fmcssn(mz01,mz01(kptimu),mz01)

!             Find SINH(2*REAL(MA)) and COSH(2*REAL(MA)).

      CALL fmadd(ma,ma,m06)
      CALL fmchsh(m06,m05,m06)

!             TANH(MA) =  SINH(2*REAL(MA)) /
!                         (COS(2*IMAG(MA))+COSH(2*REAL(MA)) +
!                         SIN(2*IMAG(MA)) /
!                         (COS(2*IMAG(MA))+COSH(2*REAL(MA)) i

      CALL fmadd(mz01(kptimu),m05,m05)
      IF (m05(2)==0) THEN
        mz01(0) = 0
        ngoal = int(real(ndsave)*alogm2) + 7
        GO TO 30
      ELSE IF (m05(1)==mexpov) THEN
        CALL fmdiv(mz01,m05,mz01(kptimu))
        CALL fmim(1,mz01)
        IF (m06(2)<0) mz01(2) = -mz01(2)
      ELSE
        CALL fmdivd(mz01,m06,m05,mz01(kptimu),mz01)
      END IF

!             Check for too much cancellation.

20    IF (ncall<=1) THEN
        ngoal = int(real(ndsave)*alogm2) + 7
      ELSE
        ngoal = int(-mxexp2)
      END IF
30    IF (mz01(0)<=ngoal .OR. mz01(kptimu)<=ngoal) THEN
        iextra = int(real(max(ngoal-mz01(0),ngoal-mz01(kptimu)))/alogm2+23.03/ &
          alogmb) + 1
        ndig = ndig + iextra
        IF (ndig>ndg2mx) THEN
          kflag = -9
          CALL zmwarn
          mz01(1) = munkno
          mz01(2) = 1
          mz01(kptimu+1) = munkno
          mz01(kptimu+2) = 1
          DO 40 j = 2, ndsave
            mz01(j+1) = 0
            mz01(kptimu+j+1) = 0
40        CONTINUE
          ndig = ndig - iextra
          mz01(0) = nint(ndig*alogm2)
          mz01(kptimu) = nint(ndig*alogm2)
          GO TO 50
        END IF
        GO TO 10
      END IF

50    maccmb = mz01(0)
      ma(0) = marz
      mz01(0) = min(maccmb,marz,maiz)
      maccmb = mz01(kptimu)
      ma(kptimu) = maiz
      mz01(kptimu) = min(maccmb,marz,maiz)
      CALL zmexit(mz01,mb,ndsave,mxsave,kasave,kovun,0)
      krad = krsave
      RETURN
    END SUBROUTINE zmtanh
    SUBROUTINE zmunpk(mp,ma)

!  MP is unpacked and the value returned in MA.

      IMPLICIT NONE

! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz), mp(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL fmunpk
! ..
      CALL fmunpk(mp,ma)
      CALL fmunpk(mp(kptimp),ma(kptimu))
      RETURN
    END SUBROUTINE zmunpk
    SUBROUTINE zmwarn

!  Called by one of the ZM routines to print a warning message
!  if any error condition arises in that routine.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC int
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
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
      ELSE IF (kflag==-8 .AND. name=='ZMOUT ') THEN
        WRITE (kw,90090)
      ELSE IF (kflag==-8 .AND. name=='ZMREAD') THEN
        WRITE (kw,90100)
      ELSE IF (kflag==-9) THEN
        WRITE (kw,90110)
        WRITE (kw,90120) ndig, ndg2mx
        WRITE (kw,90050)
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
    END SUBROUTINE zmwarn
    SUBROUTINE zmwrit(kwrite,ma)

!  Write MA on unit KWRITE under the current format.  Multi-line numbers
!  will have '&' as the last nonblank character on all but the last
!  line of the real part and the imaginary part.
!  These numbers can then be read easily using FMREAD.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC int, log10, max, min, mod, real
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kwrite
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz)
! ..
! .. Local Scalars ..
      INTEGER :: j, k, ksave, l, last, last1, last2, lb, nd, nexp
! ..
! .. External Subroutines ..
      EXTERNAL zmout
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
      CHARACTER (1) :: cmbuff(lmbuff), cmbufz(lmbufz)
      CHARACTER (6) :: namest(0:50)
! ..
! .. Common Blocks ..
      COMMON /fm/mwa, ncall, kaccsw, mxexp, mxexp2, mexpun, mexpov, munkno, &
        iunkno, runkno, mxbase, ndg2mx, spmax, dpmax, maxint, intmax, ksub
      COMMON /fmbuff/cmbuff, namest, cmchar
      COMMON /fmuser/mbase, ndig, jform1, jform2, krad, kw, ntrace, lvltrc, &
        kflag, kwarn, kround, kswide, keswch, kdebug
      COMMON /zmbuff/cmbufz
! ..
      ncall = ncall + 1
      namest(ncall) = 'ZMWRIT'
      ksave = kflag
      nd = int(real(ndig)*log10(real(mbase))) + 1
      IF (nd<2) nd = 2
      nexp = int(2.0*log10(real(mxbase))) + 6
      lb = 2*max(jform2+nexp,nd+nexp) + 3
      lb = min(lb,lmbufz)
      CALL zmout(ma,cmbufz,lb,last1,last2)
      kflag = ksave
      last = last2 + 1
      DO 10 j = 1, last2
        IF (cmbufz(last-j)/=' ' .OR. j==last2) THEN
          l = last - j
          IF (mod(l,73)/=0) THEN
            WRITE (kwrite,90000) (cmbufz(k),k=1,l)
          ELSE
            WRITE (kwrite,90000) (cmbufz(k),k=1,l-73)
            WRITE (kwrite,90010) (cmbufz(k),k=l-72,l)
          END IF
          ncall = ncall - 1
          RETURN
        END IF
10    CONTINUE
      ncall = ncall - 1
      RETURN
90000 FORMAT (4X,73A1,' &')
90010 FORMAT (4X,73A1)
    END SUBROUTINE zmwrit
    SUBROUTINE zmz2m(zval,ma)

!  MA = ZVAL

!  ZVAL is complex and is converted to ZM form.

      IMPLICIT NONE

! .. Intrinsic Functions ..
      INTRINSIC aimag, dble
! ..
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      COMPLEX :: zval
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lunpkz)
! ..
! .. Local Scalars ..
      REAL (KIND(0.0D0)) :: dz
! ..
! .. External Subroutines ..
      EXTERNAL fmdp2m, zmntr, zmntrz
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
      namest(ncall) = 'ZMZ2M '
      IF (ntrace/=0) CALL zmntrz(2,zval,1)

      dz = dble(zval)
      CALL fmdp2m(dz,ma)
      dz = dble(aimag(zval))
      CALL fmdp2m(dz,ma(kptimu))

      IF (ntrace/=0) CALL zmntr(1,ma,ma,1)
      ncall = ncall - 1
      RETURN
    END SUBROUTINE zmz2m

!  Here are the routines which work with packed ZM numbers.  All names
!  are the same as unpacked versions with 'ZM' replaced by 'ZP'.

!  To convert a program using the ZM package from unpacked calls to
!  packed calls make these changes to the program:
!  '(0:LUNPKZ)' to '(0:LUNPKZ)' in dimensions.
!  'CALL ZM' to 'CALL ZP'

    SUBROUTINE zpabs(ma,mbfm)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mbfm(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack, zmabs, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mpa(0:lunpck), mpb(0:lunpck), mpc(0:lunpck), mx(0:lunpkz), &
        my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mpa, mpb, mpc
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmabs(mx,mpa)
      CALL fmpack(mpa,mbfm)
      RETURN
    END SUBROUTINE zpabs
    SUBROUTINE zpacos(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmacos, zmpack, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmacos(mx,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zpacos
    SUBROUTINE zpadd(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz), mc(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmadd, zmpack, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmunpk(mb,my)
      CALL zmadd(mx,my,mx)
      CALL zmpack(mx,mc)
      RETURN
    END SUBROUTINE zpadd
    SUBROUTINE zpaddi(ma,integ)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: integ
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmaddi, zmpack, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmaddi(mx,integ)
      CALL zmpack(mx,ma)
      RETURN
    END SUBROUTINE zpaddi
    SUBROUTINE zparg(ma,mbfm)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mbfm(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack, zmarg, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mpa(0:lunpck), mpb(0:lunpck), mpc(0:lunpck), mx(0:lunpkz), &
        my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mpa, mpb, mpc
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmarg(mx,mpa)
      CALL fmpack(mpa,mbfm)
      RETURN
    END SUBROUTINE zparg
    SUBROUTINE zpasin(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmasin, zmpack, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmasin(mx,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zpasin
    SUBROUTINE zpatan(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmatan, zmpack, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmatan(mx,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zpatan
    SUBROUTINE zpchsh(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz), mc(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmchsh, zmpack, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmchsh(mx,mx,my)
      CALL zmpack(mx,mb)
      CALL zmpack(my,mc)
      RETURN
    END SUBROUTINE zpchsh
    SUBROUTINE zpcmpx(mafm,mbfm,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: mafm(0:lpack), mbfm(0:lpack), mc(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL fmunpk, zmcmpx, zmpack
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mpa(0:lunpck), mpb(0:lunpck), mpc(0:lunpck), mx(0:lunpkz), &
        my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mpa, mpb, mpc
      COMMON /zmpck/mx, my
! ..

      CALL fmunpk(mafm,mpa)
      CALL fmunpk(mbfm,mpb)
      CALL zmcmpx(mpa,mpb,mx)
      CALL zmpack(mx,mc)
      RETURN
    END SUBROUTINE zpcmpx
    SUBROUTINE zpconj(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmconj, zmpack, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmconj(mx,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zpconj
    SUBROUTINE zpcos(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmcos, zmpack, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmcos(mx,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zpcos
    SUBROUTINE zpcosh(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmcosh, zmpack, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmcosh(mx,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zpcosh
    SUBROUTINE zpcssn(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz), mc(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmcssn, zmpack, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmcssn(mx,mx,my)
      CALL zmpack(mx,mb)
      CALL zmpack(my,mc)
      RETURN
    END SUBROUTINE zpcssn
    SUBROUTINE zpdiv(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz), mc(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmdiv, zmpack, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmunpk(mb,my)
      CALL zmdiv(mx,my,mx)
      CALL zmpack(mx,mc)
      RETURN
    END SUBROUTINE zpdiv
    SUBROUTINE zpdivi(ma,integ,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: integ
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmdivi, zmpack, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmdivi(mx,integ,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zpdivi
    SUBROUTINE zpeq(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL fpeq
! ..

      CALL fpeq(ma,mb)
      CALL fpeq(ma(kptimp),mb(kptimp))
      RETURN
    END SUBROUTINE zpeq
    SUBROUTINE zpequ(ma,mb,nda,ndb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: nda, ndb
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL fpequ
! ..

      CALL fpequ(ma,mb,nda,ndb)
      CALL fpequ(ma(kptimp),mb(kptimp),nda,ndb)
      RETURN
    END SUBROUTINE zpequ
    SUBROUTINE zpexp(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmexp, zmpack, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmexp(mx,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zpexp
    SUBROUTINE zpform(form1,form2,ma,string)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (*) :: form1, form2, string
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmform, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmform(form1,form2,mx,string)
      RETURN
    END SUBROUTINE zpform
    SUBROUTINE zpfprt(form1,form2,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (*) :: form1, form2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmfprt, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmfprt(form1,form2,mx)
      RETURN
    END SUBROUTINE zpfprt
    SUBROUTINE zp2i2m(integ1,integ2,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: integ1, integ2
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zm2i2m, zmpack
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zm2i2m(integ1,integ2,mx)
      CALL zmpack(mx,ma)
      RETURN
    END SUBROUTINE zp2i2m
    SUBROUTINE zpi2m(integ,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: integ
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmi2m, zmpack
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmi2m(integ,mx)
      CALL zmpack(mx,ma)
      RETURN
    END SUBROUTINE zpi2m
    SUBROUTINE zpimag(ma,mbfm)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mbfm(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack, zmimag, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mpa(0:lunpck), mpb(0:lunpck), mpc(0:lunpck), mx(0:lunpkz), &
        my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mpa, mpb, mpc
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmimag(mx,mpa)
      CALL fmpack(mpa,mbfm)
      RETURN
    END SUBROUTINE zpimag
    SUBROUTINE zpint(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmint, zmpack, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmint(mx,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zpint
    SUBROUTINE zpinp(line,ma,la,lb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: la, lb
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz)
      CHARACTER (1) :: line(lb)
! ..
! .. External Subroutines ..
      EXTERNAL zminp, zmpack
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zminp(line,mx,la,lb)
      CALL zmpack(mx,ma)
      RETURN
    END SUBROUTINE zpinp
    SUBROUTINE zpipwr(ma,integ,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: integ
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmipwr, zmpack, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmipwr(mx,integ,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zpipwr
    SUBROUTINE zplg10(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmlg10, zmpack, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmlg10(mx,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zplg10
    SUBROUTINE zpln(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmln, zmpack, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmln(mx,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zpln
    SUBROUTINE zpm2i(ma,integ)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: integ
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmm2i, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmm2i(mx,integ)
      RETURN
    END SUBROUTINE zpm2i
    SUBROUTINE zpm2z(ma,zval)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      COMPLEX :: zval
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmm2z, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmm2z(mx,zval)
      RETURN
    END SUBROUTINE zpm2z
    SUBROUTINE zpmpy(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz), mc(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmmpy, zmpack, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmunpk(mb,my)
      CALL zmmpy(mx,my,mx)
      CALL zmpack(mx,mc)
      RETURN
    END SUBROUTINE zpmpy
    SUBROUTINE zpmpyi(ma,integ,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: integ
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmmpyi, zmpack, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmmpyi(mx,integ,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zpmpyi
    SUBROUTINE zpnint(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmnint, zmpack, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmnint(mx,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zpnint
    SUBROUTINE zpout(ma,line,lb,last1,last2)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: last1, last2, lb
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz)
      CHARACTER (1) :: line(lb)
! ..
! .. External Subroutines ..
      EXTERNAL zmout, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmout(mx,line,lb,last1,last2)
      RETURN
    END SUBROUTINE zpout
    SUBROUTINE zpprnt(ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmprnt, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmprnt(mx)
      RETURN
    END SUBROUTINE zpprnt
    SUBROUTINE zppwr(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz), mc(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmpack, zmpwr, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmunpk(mb,my)
      CALL zmpwr(mx,my,mx)
      CALL zmpack(mx,mc)
      RETURN
    END SUBROUTINE zppwr
    SUBROUTINE zpread(kread,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kread
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmpack, zmread
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmread(kread,mx)
      CALL zmpack(mx,ma)
      RETURN
    END SUBROUTINE zpread
    SUBROUTINE zpreal(ma,mbfm)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mbfm(0:lpack)
! ..
! .. External Subroutines ..
      EXTERNAL fmpack, zmreal, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mpa(0:lunpck), mpb(0:lunpck), mpc(0:lunpck), mx(0:lunpkz), &
        my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /fmpck/mpa, mpb, mpc
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmreal(mx,mpa)
      CALL fmpack(mpa,mbfm)
      RETURN
    END SUBROUTINE zpreal
    SUBROUTINE zprpwr(ma,ival,jval,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: ival, jval
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmpack, zmrpwr, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmrpwr(mx,ival,jval,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zprpwr
    SUBROUTINE zpset(nprec)
! .. Scalar Arguments ..
      INTEGER :: nprec
! ..
! .. External Subroutines ..
      EXTERNAL zmset
! ..

      CALL zmset(nprec)
      RETURN
    END SUBROUTINE zpset
    SUBROUTINE zpsin(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmpack, zmsin, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmsin(mx,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zpsin
    SUBROUTINE zpsinh(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmpack, zmsinh, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmsinh(mx,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zpsinh
    SUBROUTINE zpsqr(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmpack, zmsqr, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmsqr(mx,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zpsqr
    SUBROUTINE zpsqrt(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmpack, zmsqrt, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmsqrt(mx,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zpsqrt
    SUBROUTINE zpst2m(string,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      CHARACTER (*) :: string
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmpack, zmst2m
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmst2m(string,mx)
      CALL zmpack(mx,ma)
      RETURN
    END SUBROUTINE zpst2m
    SUBROUTINE zpsub(ma,mb,mc)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz), mc(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmpack, zmsub, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmunpk(mb,my)
      CALL zmsub(mx,my,mx)
      CALL zmpack(mx,mc)
      RETURN
    END SUBROUTINE zpsub
    SUBROUTINE zptan(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmpack, zmtan, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmtan(mx,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zptan
    SUBROUTINE zptanh(ma,mb)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz), mb(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmpack, zmtanh, zmunpk
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmtanh(mx,mx)
      CALL zmpack(mx,mb)
      RETURN
    END SUBROUTINE zptanh
    SUBROUTINE zpwrit(kwrite,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      INTEGER :: kwrite
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmunpk, zmwrit
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmunpk(ma,mx)
      CALL zmwrit(kwrite,mx)
      RETURN
    END SUBROUTINE zpwrit
    SUBROUTINE zpz2m(zval,ma)
      IMPLICIT NONE
! .. Parameters ..
      INTEGER, PARAMETER :: nbits = 64, ndigmx = 256
      INTEGER, PARAMETER :: lpack = (ndigmx+1)/2 + 1
      INTEGER, PARAMETER :: lpackz = 2*lpack + 1
      INTEGER, PARAMETER :: lunpck = (6*ndigmx)/5 + 20
      INTEGER, PARAMETER :: lunpkz = 2*lunpck + 1
      INTEGER, PARAMETER :: kptimp = lpack + 1
      INTEGER, PARAMETER :: kptimu = lunpck + 1
      INTEGER, PARAMETER :: ljsums = 8*(lunpck+2)
      INTEGER, PARAMETER :: lmbuff = ((lunpck+3)*(nbits-1)*301)/2000 + 6
      INTEGER, PARAMETER :: lmbufz = 2*lmbuff + 10
      INTEGER, PARAMETER :: lmwa = 2*lunpck
! ..
! .. Scalar Arguments ..
      COMPLEX :: zval
! ..
! .. Array Arguments ..
      REAL (KIND(0.0D0)) :: ma(0:lpackz)
! ..
! .. External Subroutines ..
      EXTERNAL zmpack, zmz2m
! ..
! .. Arrays in Common ..
      REAL (KIND(0.0D0)) :: mx(0:lunpkz), my(0:lunpkz)
! ..
! .. Common Blocks ..
      COMMON /zmpck/mx, my
! ..

      CALL zmz2m(zval,mx)
      CALL zmpack(mx,ma)
      RETURN
!             End of the ZM package.
    END SUBROUTINE zpz2m
